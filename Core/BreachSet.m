classdef BreachSet < BreachStatus
    % BreachSet Defines an API to manipulate parameters and traces.
    %    This is the base class for BreachSystem objects. BreachSet
    %    instances are usually also BreachSystem or some derived class
    %    instances, so the purpose of this class is mainly to separate
    %    methods for defining and manipulate parameters and traces data.
    %
    % BreachSet Properties
    %   ParamRanges              - ranges of possible values for each parameter - determines the parameter sampling domain
    %   SignalRanges             - ranges of values taken by each signal variable
    %   AppendWhenSample=false   - when true, sampling appends new param vectors, otherwise replace.
    %
    % BreachSet Methods
    %   GetParam            - get values for parameters, given their names
    %   SetParam            - set values for parameters, given their names
    %   GetParamRanges      - get ranges for parameters, given their names and ranges
    %   SetParamRanges      - set ranges for parameters, given their names and ranges
    %   GetDomain
    %   SetDomain   
    %   SampleDomain 
    %   GridSample          - creates parameter vectors in ParamRanges on a regularly spaced grid
    %   CornerSample        - creates parameter vectors from the corners of ParamRanges
    %   QuasiRandomSample   - uniformly samples n parameter vectors in ParamRanges
    %   Concat              - concatenates parameter vectors and traces of another compatible BreachSet
    %   PrintParams         - display the number of parameter vectors, parameter names and ranges
    %   PrintSignals        - display the number of traces, names of signal variables and ranges
    %   PlotParams          - plots parameter vectors
    %   PlotSignals         - plots signals vs time
    %   PlotSigPortrait     - plots signals portrait
    
    properties
        P % legacy parameter structure - contains points data (in P.pts) and traces (P.traj) and many other fields whose purpose is slowly falling into oblivion
        Domains = BreachDomain('double', [])
        % will replace SignalRanges eventually (?)
        SignalRanges % ranges of values taken by each signal variable      
        AppendWhenSample=false % when true, sampling appends new param vectors, otherwise replace.
    end
    
    methods (Hidden=true)
        function P = GetP(this)
            % Get the legacy parameter set structure
            P = this.P;
        end
        
        function SetP(this, P)
            % Get the legacy parameter set structure
            this.P = P;
            this.CheckinDomain();
        end
    end
    
    methods
        %% Constructor
        function this = BreachSet(Sys, params, ranges)
            % BreachSet constructor from a legacy P or Sys, parameter names and ranges
            
            switch nargin
                case 0
                    return;
                case 1
                    this.P = CreateParamSet(Sys);
                case 2
                    this.P = CreateParamSet(Sys, params);
                case 3
                    this.P = CreateParamSet(Sys, params, ranges);
            end
            this.CheckinDomain();
            
        end
        
        %% Domains
        function SetDomain(this, params, type, domain)
            % SetDomain
            
            if nargin<4
                domain =ones(numel(params),0);
            end
            
            if ischar(params)
                params = {params};
            end
            
            if isa(type,'BreachDomain')
                if numel(type)==1
                    type = repmat(type,1, numel(params));
                end
                if numel(type) ~= numel(params)
                    error('number of types and parameters or signals mismatch.')
                end
            elseif ischar(type)
                com_type = type; 
                type = repmat(BreachDomain(), 1, numel(params));
                for ip = 1:numel(params)
                    type(ip) = BreachDomain(com_type, domain(ip,:));  
                end
            elseif iscell(type)
                cell_type = type;
                type = repmat(BreachDomain(), 1, numel(params));
                for ip = 1:numel(params)
                    if isa(cell_type{ip}, 'BreachDomain')
                        type(ip) = type{ip};
                    else 
                        type(ip) = BreachDomain(cell_type{ip}, domain(ip,:));
                    end
               end
             end
             
             for ip = 1:numel(params)
                 param = params{ip};
                 [idx, found] = FindParam(this.P, param);
                 if any(found==0)
                     error('Parameter or signal not found.');
                 end
                 
                 if isa(type(ip),'BreachDomain')
                     this.Domains(idx) = type(ip);
                 else
                     this.Domains(idx) = BreachDomain(type{ip,:}, domain(ip,:));
                     if ~isempty(domain(ip,:))&&idx>this.P.DimX&&~isequal(type(ip,:), 'enum')
                         this.SetParamRanges(idx, [this.Domains(idx).domain(1),this.Domains(idx).domain(2)]);
                     end
                 end
             end
               this.CheckinDomain();
        end
        
        function dom = GetDomain(this, param)
            % GetDomain
            idx = FindParam(this.P, param);
            for i=1:numel(idx)
                if numel(this.Domains)<idx(i)
                    dom(i)= BreachDomain();
                else
                    dom(i) = this.Domains(idx(i));
                end
            end
            if ~exist('dom','var')
                dom = [];
            end
        end
        
        function CheckinDomain(this)
            % CheckinDomain
            this.CheckinDomainParam();
            this.CheckinDomainTraj();
        end
        
        function CheckinDomainParam(this)
            pts = this.P.pts;
            if numel(this.Domains)< size(pts,1)
                this.Domains(size(pts,1)) = BreachDomain();
            end
            for  i =this.P.DimX+1:size(pts,1)
                if ~isequal(this.Domains(i).type, 'double')||~isempty(this.Domains(i).domain)
                    pts(i,:) = this.Domains(i).checkin(pts(i,:));
                end
            end
            
            this.P.pts = pts;
            
        end
        
        function CheckinDomainTraj(this)
            % CheckinDomainTraj
            % checks trajectories
            if this.hasTraj()
                for itraj = 1:numel(this.P.traj)
                    for  i=1:this.P.DimX
                        if ~isempty(this.Domains(i).domain)
                            this.P.traj{itraj}.X(i,:) = this.Domains(i).checkin(this.P.traj{itraj}.X(i,:));
                        end
                    end
                end
            end
        end
        
        %%  Param 
        function SetParam(this, params, values, is_spec_param)
        % BreachSet.SetParam(params, values,  is_spec_param) sets values to
        % parameters listed in params. If the set contains only one sample,
        % creates as many sample as there are values. If the set has
        % several samples and there is only one value, set this value to
        % all samples. Otherwise, returns an error. 
            
            ip = FindParam(this.P, params);
            i_not_sys = find(ip>this.P.DimP);
            if ~isempty(i_not_sys)
                if iscell(params)
                    nparam = params{i_not_sys(1)};
                else
                    nparam = params;
                end
                if ~exist('is_spec_param','var')||isequal(is_spec_param,false)
                    warning('SetParam:param_not_in_list',['Parameter ' nparam ' was set but is not a system parameter.' ...
                        ' If this is intended to be a spec. parameter, consider using SetParamSpec instead.']);
                end
            end
            
            num_params = numel(ip);
            
            if size(values, 1)~= num_params
                if size(values,1) == 1 && size(values,2 ) == num_params
                    values = values';
                else
                    error('SetParam:wrong_arguments_size', 'Dimension mismatch between values and parameters.');
                end
            end
            
            num_pts  =  size(this.P.pts,2);
            num_values = size(values, 2); 
            
            if num_values==1 || num_values == num_pts
                this.P = SetParam(this.P, params, values);
            elseif num_pts==1    % note in this case, we have to remove traces ( or see if maybe not, ) 
                this.P = Sselect(SPurge(this.P),1);
                this.P.pts = repmat(this.P.pts,1, size(values, 2));
                this.P.epsi= repmat(this.P.epsi,1, size(values, 2));
                this.P.selected = zeros(1, size(values, 2));
                this.P = SetParam(this.P, ip, values);
            else
                error('SetParam:wrong_arguments_size', 'Dimension mismatch between values and parameters.');
            end
        end
                
        function SetParamSpec(this, params, values, ignore_sys_param)
            ip = FindParam(this.P, params);
            if all(ip>this.P.DimP)
                this.P = SetParam(this.P, params, values);
                % adds Domains 
                for i = ip
                    if size(this.Domains, 2)<i
                        this.Domains(i) = BreachDomain();
                    end
                end
            elseif ~exist('ignore_sys_param', 'var')||ignore_sys_param==false 
                error('Attempt to modify a system parameter - use SetParam instead.');
            end
        end
        
        function values = GetParam(this, params, ip)
            values = GetParam(this.P,params);
            if exist('ip', 'var')
                values = values(:, ip);
            end
        end
        
        function ResetParamSet(this)
            % ResetParamSet remove samples and keeps one in the domain  
        
            this.P = SPurge(this.P);
            % find non empty domains
            ipr = cellfun(@(c)(~isempty(c)), {this.Domains.domain});            
            if (isempty(find(ipr,1)))
                this.P = CreateParamSet(this.P);
                this.P.epsi(:,:)=0;
            else 
                ranges = cat( 1 , this.Domains.domain );
                this.P = CreateParamSet( this.P , this.P.ParamList(ipr) , ranges);
                this.CheckinDomain();
            end
        end
        
        
        %% Get and Set param ranges
        function SetParamRanges(this, params, ranges)
        % SetParamRanges set intervals for parameters (set domains as
        % bounded 'double' or 'int' if it is already an 'int') 
            i_not_found= [];
            if ~isnumeric(params)
                [i_params, res]= FindParam(this.P, params);
                i_not_found = find(res==0);
            else
                i_params=params;
            end
            
            if size(ranges,1) == 1
                ranges = repmat(ranges, numel(i_params),1);
            end
            
            if ~isempty(i_not_found)
                if iscell(params)
                    param_not_found = params{i_not_found(1)};
                else
                    param_not_found = params;
                end
                error('SetParamRanges:param_not_found', ['Parameter ' param_not_found ' not found.']);
            end
            
            % Set domain 
            for ip = 1:numel(i_params)
                   type = this.Domains(i_params(ip)).type;
                   if isequal(type, 'enum')||isequal(type,'bool')
                      warning('SetParamRanges:enum_or_bool', 'Use SetDomain for enum or bool types.' );
                   else
                       this.Domains(i_params(ip)) = BreachDomain(type, ranges(ip,:));
                   end
            end
            
            this.P.dim = i_params;
            
        end
        
        function ranges = GetParamRanges(this, params)
        % GetParamRanges 
            i_params = FindParam(this.P, params);
            ranges= zeros(numel(params),2);
            ranges(:,1) = -inf;
            ranges(:,2) = inf;

            for ip = 1:numel(i_params)
                if (~isempty(this.Domains(i_params(ip)).domain))
                    ranges(ip,:) =  this.Domains(i_params(ip)).domain;
                end
            end
            
        end
        
        function ResetEpsi(this)
            % Set Param ranges around individual parameter vectors to zero
            this.P.epsi(:,:) = 0;
        end
                
        function traces = GetTraces(this)
            % Get computed trajectories
            traces= [];
            if isfield(this.P,'traj')
                traces = this.P.traj;
            end
        end
        
        function val = UpdateSignalRanges(this)
            % Update ranges for variables from trajectories in P
            if isfield(this.P, 'traj')
                if isempty(this.SignalRanges)
                    this.SignalRanges = ones(this.Sys.DimX,2);
                    minX = +inf*ones(this.Sys.DimX,1);
                    maxX = -inf*ones(this.Sys.DimX,1);
                else
                    minX = this.SignalRanges(:,1);
                    maxX = this.SignalRanges(:,2);
                end
                val=inf;
                for itraj = 1:numel(this.P.traj)
                    traj = this.P.traj{itraj};
                    traj_maxX = max(traj.X,[], 2);
                    traj_minX = min(traj.X,[], 2);
                    dist_maxX = min(maxX-traj_maxX);
                    dist_minX = min(traj_minX-minX);
                    val= min( [val dist_maxX dist_minX] );
                    minX = min([traj_minX minX],[],2);
                    maxX = max([traj_maxX maxX],[],2);
                end
                this.SignalRanges = [minX, maxX];
                this.P.SignalRanges = this.SignalRanges; % duplicate - never good for sure...
            end
            
        end
        
        function SigNames = GetSignalNames(this)
            % Get signal names - same as GetSignalList
            SigNames = this.P.ParamList(1:this.P.DimX);
        end
        
        function signals = GetSignalList(this)
            % GetSignalList returns signal names
            signals = this.P.ParamList(1:this.P.DimX);
        end
        
        function params = GetParamList(this)
            % GetParamList returns parameter names
            params = this.P.ParamList(this.P.DimX+1:end);
        end
        
        function sys_params = GetSysParamList(this)
            % GetSysParamList returns system parameter names
            sys_params = this.P.ParamList(this.P.DimX+1:this.P.DimP);
        end
        
        function prop_params = GetPropParamList(this)
            % GetSysParamList returns system parameter names
            prop_params = this.P.ParamList(this.P.DimP+1:end);
        end
        
        
        function X = GetSignalValues(this, iX, t)
            % Get signal values - in case of several trajectories, return cell array
            if (~isfield(this.P,'traj'))
                error('GetTrajValues:NoTrajField','Compute/import trajectories first.')
            end
            
            if ischar(iX) || iscell(iX)
                iX = FindParam(this.P, iX);
            end
            
            nb_traj = numel(this.P.traj);
            X = cell(nb_traj,1);
            for i_traj = 1:nb_traj
                if (~exist('t','var'))
                    X{i_traj} = this.P.traj{i_traj}.X(iX,:);
                else
                    X{i_traj} = interp1(this.P.traj{i_traj}.time, this.P.traj{i_traj}.X(iX,:)',t)';
                    if numel(iX)==1
                        X{i_traj} = X{i_traj}';
                    end
                end
                
            end
            if nb_traj==1
                X = X{1};
            end
        end
        
        function h = PlotSignals(this, varargin)
            % Plot signals
            if (~isfield(this.P,'traj'))
                error('No signal to plot. Use Sim command first.')
            end
            
            gca;
            h = SplotVar(this.P, varargin{:});
            
        end
        
        function h = PlotSigPortrait(this, varargin)
            % Plot signals, phase portrait
            if (~isfield(this.P,'traj'))
                error('No signal to plot. Use Sim command first.')
            end
            
            gca;
            SplotTraj(this.P, varargin{:});
        end
        
        %% Sampling
        function SampleDomain(this, params, num_samples, method, opt_multi)
            % BreachSet.SampleDomain generic sampling function 
            %
            % B.SampleDomain('p', 5) creates 5 samples drawn randomly
            % from the domain of parameter p. Do nothing if domain is
            % empty (which it is by default).
            %
            % B.SampleDomain('p','all') create samples enumerating all
            % possible values in the domain of 'p'.
            %
            % B.SampleDomain('p',5,'grid') draw samples regularly spaced, including corners
            %
            % B.SampleDomain({'p1','p2'},5,'grid') 5x5 grid
            %
            % B.SampleDomain(...,...,..., 'replace') default
            %
            % B.SampleDomain(...,...,..., 'append')
            %
            % B.SampleDomain(...,...,..., 'combine')
            %
            
            % process parameters
            if ischar(params)
                params = {params};
            end
            
            idx_param = FindParam(this.P, params);
            
            % if we have traces, we'll need to save and restore them 
            saved_traj = false; 
            if this.hasTraj()
                saved_traj = true;
                P0 = this.P;
            end
            
            domains = this.Domains(idx_param);
            domains =  num2cell(domains); % convert array to cell

            if ~exist('method')||isempty(method)
                method = 'rand';
            end
            
            if ~exist('opt_multi')
                opt_multi='replace'; %
            end
            
            if isequal(method, 'quasi-random')&&(~isequal(num_samples, 'all')|| (iscell(num_samples)&&any(strcmp(num_samples, 'all'))))
            % FIXME can do better
                Pold = this.P;
                Appold = this.AppendWhenSample;
                this.P.dim = idx_param;
                this.AppendWhenSample =0;
                this.QuasiRandomSample(prod(num_samples));
                x = this.GetParam(idx_param);
                this.P = Pold;
            else
                x = sample(domains{:}, num_samples, method);
            end
          
            % combine (or not) with others
            num_old = this.GetNbParamVectors();
            switch opt_multi
                case 'replace'
                    % if there are multiple samples already, discard them
                    % and use default p
                    if num_old==1
                        this.SetParam(params, x, true)
                    else
                        this.ResetParamSet();
                        this.SetParam(params, x, true);
                    end
                    
                case 'append'  % FIXME deep copy here maybe not smartest   
                    Btmp = this.copy();
                    if num_old==1
                        Btmp.SetParam(params, x, true)
                    else
                        Btmp.ResetParamSet();
                        Btmp.SetParam(params, x, true);
                        this.Concat(Btmp);
                    end
                    
                case 'combine'
                    num_new = numel(x);
                    idx = N2Nn(2, [num_old num_new]);
                    pts = this.P.pts(:, idx(1,:));
                    pts(idx_param,:) = x(:, idx(2,:));
                    this.ResetParamSet();
                    this.SetParam(1:size(pts,1),pts, true);
            end
            
            % restore traj if needed 
            if saved_traj
                this.P = Pfix_traj_ref(this.P, P0);
            end
            
        end
        
    
        %% Legacy sampling
        function GridSample(this, delta)
            % BreachSet.GridSample(num_samples) sample all bounded domain
            % with a grid of num_samples elements. num_samples may be a
            % scalar or an array of same dimension as the number of bounded
            % domains.
            bnd_params = this.GetBoundedDomains();
            if this.AppendWhenSample
               this.SampleDomain(bnd_params, delta, 'grid', 'append');
            else
               this.SampleDomain(bnd_params, delta, 'grid');
            end
            
        end
        
        % Get corners of parameter domain
        function CornerSample(this)
            if this.AppendWhenSample
                Pold = this.P;
            end
            
            this.ResetParamSet();
            newP = this.P;
            newP.epsi = 2*newP.epsi;
            newP = Refine(newP,2);
            newP.epsi = newP.epsi/2;
            
            if this.AppendWhenSample
                this.P = SConcat(Pold, newP);
            else
                this.P = newP;
            end
            this.CheckinDomainParam();
            
        end
        
        function QuasiRandomSample(this, nb_sample, step)
            % Quasi-Random Sampling
                     
            if this.AppendWhenSample
                Pold = this.P;
            end
            
            this.ResetParamSet();
            if nargin==3
                newP = QuasiRefine(this.P,nb_sample, step);
            else
                newP = QuasiRefine(this.P, nb_sample);
            end
            
            if this.AppendWhenSample
                this.P = SConcat(Pold, newP);
            else
                this.P = newP;
            end
            this.CheckinDomainParam();
        end
        
        % Get the number of param vectors - -1 means P is empty
        function nb_pts = GetNbParamVectors(this)
            if isempty(this.P)
                nb_pts = -1;
            else
                nb_pts= size(this.P.pts,2);
            end
        end
        
        %% Concatenation - needs some additional compatibility checks...
        function Concat(this, other)
            this.P = SConcat(this.P, other.P);
        end
        
        %% Plot parameters
        function PlotParams(this, varargin)
            % Plot parameters
            gca;
            P = DiscrimPropValues(this.P);
            SplotPts(P, varargin{:});
        end
        
        %% Coverage
        function [cnt, grd1, grd2] = GetSignalCoverage(this,sigs, delta1,delta2)
            % 1d or 2d
            X = this.GetSignalValues(sigs);
            
            switch (numel(sigs))
                case 1
                    [cnt, grd1] = cover(X,delta1);
                case 2
                    [cnt, grd1, grd2] = cover2d(X,delta1,delta2);
                otherwise
                    error('Coverage for more than 2 signals is not supported');
            end
        end
        
        function [cnt, grd1, grd2] = PlotSignalCoverage(this,sigs, delta1,delta2)
            % 1d or 2d
            X = this.GetSignalValues(sigs);
            switch (numel(sigs))
                case 1
                    [cnt, grd1] = cover(X,delta1);
                    X = grd1+delta1/2; % centers of bins
                    bar(X,cnt)
                    set(gca, 'XTick', round([grd1 grd1(end)+delta1]))
                    title('number of samples per bin')
                    xlabel(sigs{1});
                case 2
                    [cnt, grd1, grd2] = cover2d(X,delta1,delta2);
                    X = grd1+delta1/2; % centers of bins
                    Y = grd2+delta2/2; % centers of bins
                    figure;
                    b= bar3(cnt');
                    set(gca,'XTickLabel', round([X X(end)+delta1]) )
                    set(gca,'YTickLabel', round([Y Y(end)+delta2]) )
                    %surface(grd2, grd1, cnt);
                    xlabel(sigs{1});
                    ylabel(sigs{2});
                    colormap(jet);
                    colorbar;
                    
                    for k = 1:length(b)
                        zdata = b(k).ZData;
                        b(k).CData = zdata;
                        b(k).FaceColor = 'interp';
                    end
                    title('num. samples per grid element')
                otherwise
                    error('Coverage for more than 2 signals is not supported');
            end
        end
        
        function PlotDomain(this, params)
        % PlotDomain 
        
        % default style
        col = [0 0 1];
        alpha = 0.03;
        pts_style = 'sb';
       
        
        % default params
         if ~exist('params', 'var') || isempty('params')
              params =this.GetBoundedDomains(); 
         end
         
         if ischar(params)
             params = {params};
         end
         
         switch numel(params)
             case 1 % one domain
                  
                domain1 = this.GetDomain(params);
                plotxdomain(domain1,0); 
                %  Labels 
                 xlabel(params{1}, 'Interpreter', 'none');
                set(gca, 'YTick', []);
                grid on;
                xlim = get(gca, 'XLim');
                dxlim = xlim(2) - xlim(1);
                set(gca, 'XLim', [xlim(1)-dxlim/10, xlim(2)+dxlim/10]);
                
             case 2 % two domains 
                 domain1 = this.GetDomain(params{1});
                 domain2 = this.GetDomain(params{2});
                 if isempty(domain1.enum)    % domain1 is dense   
        
                     if isempty(domain2.enum)  % domain2 is dense
                       start = [domain1.domain(1), domain2.domain(1)];
                       sz = [domain1.domain(2) - domain1.domain(1),domain2.domain(2) - domain2.domain(1)];
                       rect(start, sz, col,alpha);
                   else % domain1 dense and domain2 not dense    
                        y = domain2.sample_all();
                        hold on;
                        for i_y = 1:numel(y)
                              plotxdomain( domain1, y(i_y) );                   
                        end
                     end
                     
                 else  % domain1 not dense
                   
                     if isempty(domain2.enum)  % domain2 is dense
                        x = domain1.sample_all();
                        hold on;
                        for i_x = 1:numel(x)
                              plotydomain( domain1, x(i_x) );                   
                        end
                     else % domain1 and domain2 not dense    
                        X = sample( domain1, domain2, 'all');
                        plot(X(1,:), X(2,:), pts_style);
                     end
                 end
                 
                 xlabel(params{1}, 'Interpreter', 'none');
                 ylabel(params{2}, 'Interpreter', 'none');
                 grid on;
                xlim = get(gca, 'XLim');
                dxlim = xlim(2) - xlim(1);
                set(gca, 'XLim', [xlim(1)-dxlim/10, xlim(2)+dxlim/10]);
                ylim = get(gca, 'YLim');
                dylim = ylim(2) - ylim(1);
                set(gca, 'YLim', [ylim(1)-dylim/10, ylim(2)+dylim/10]);
                 
             case 3
                 
                 domain1 = this.GetDomain(params{1});
                 domain2 = this.GetDomain(params{2});
                 domain3 = this.GetDomain(params{3});
                 
                     if ~isempty(domain1.enum)&&~isempty(domain2.enum)&&~isempty(domain3.enum)  
                         X = sample( domain1, domain2, domain3, 'all');
                        plot3(X(1,:), X(2,:), X(3,:),  pts_style);
                 
                     else
                       start = [domain1.domain(1), domain2.domain(1), domain3.domain(1)];
                       sz = [domain1.domain(2) - domain1.domain(1),domain2.domain(2) - domain2.domain(1),domain3.domain(2) - domain3.domain(1)];
                       voxel(start, sz, col, alpha);    
                     end
                     % TODO missing cases combining enum and dense
                         
                 xlabel(params{1}, 'Interpreter', 'none');
                 ylabel(params{2}, 'Interpreter', 'none');
                 zlabel(params{3}, 'Interpreter', 'none');
                 grid on;
                xlim = get(gca, 'XLim');
                dxlim = xlim(2) - xlim(1);
                set(gca, 'XLim', [xlim(1)-dxlim/10, xlim(2)+dxlim/10]);
                ylim = get(gca, 'YLim');
                dylim = ylim(2) - ylim(1);
                set(gca, 'YLim', [ylim(1)-dylim/10, ylim(2)+dylim/10]);
                zlim = get(gca, 'ZLim');
                dzlim = zlim(2) - zlim(1);
                set(gca, 'ZLim', [zlim(1)-dzlim/10, zlim(2)+dzlim/10]);
                view([ 45 45 ]);
         end
         
            % 1d x direction
            function  plotxdomain(dom, y0)
                 if isempty(dom.enum) %  it's a dense box
                     start = [dom.domain(1), y0];
                     sz = [dom.domain(2) - dom.domain(1), 0];
                     rect(start, sz, col, alpha);
                 else
                     x = dom.sample_all(); 
                     plot(x, 0*x+y0, pts_style);
                 end                
            end
            
            % 1d y-direction
            function  plotydomain(dom, x0)
                 if isempty(dom.enum) %  it's a dense box
                     start = [x0, dom.domain(1)];
                     sz = [0, dom.domain(2) - dom.domain(1)];
                     rect(start, sz, col, alpha);
                 else
                     y = dom.sample_all(); 
                     plot(0*y+x0, y, pts_style);
                 end                
            end
        end
        
        function [ params, ipr]  = GetBoundedDomains(this)
            % GetNonEmptyDomains
            ipr = cellfun(@(c)(~isempty(c)), {this.Domains.domain});
            params =   this.P.ParamList(ipr);
        end
        
        %% Printing
        function PrintSignals(this)
            if isempty(this.SignalRanges)
                disp( 'Signals:')
                disp( '-------')
                for isig = 1:this.P.DimX
                    fprintf('%s\n', this.P.ParamList{isig});
                end
            else
                
                disp('-------')
                for isig = 1:this.P.DimX
                    fprintf('%s %s\n', this.P.ParamList{isig}, this.Domains(isig).short_disp());
                end
            end
            disp(' ')
        end
        
        function PrintParams(this)
            nb_pts= this.GetNbParamVectors();
            if (nb_pts<=1)
                disp('Parameters:')
                disp('----------')
                for ip = this.P.DimX+1:numel(this.P.ParamList)
                    fprintf('%s=%g       %s',this.P.ParamList{ip},this.P.pts(ip,1), this.Domains(ip).short_disp(1));
                    fprintf('\n');
                end
            else
                fprintf('Parameters (%d vectors):\n',nb_pts);
                disp('-------------------------');
                for ip = this.P.DimX+1:numel(this.P.ParamList)
                        fprintf('%s     %s\n',this.P.ParamList{ip}, this.Domains(ip).short_disp(1));
                end
            end
            
            disp(' ')
        end
        
        %% Misc
        
        % Warning handler
        function WarningResetP(this, fname)
            if this.GetNbParamVectors()>1
                warning('BreachSet:warning_will_reset_pts',['This set contains more than one parameter vector or traces - the function ' fname ' will likely erase them.'])
            end
        end
        
        function bl = hasTraj(this)
            % checks if this has traces
            bl=isfield(this.P, 'traj');
        end
        
        function Reset(this) 
            % Resets 
            this.P = Sselect(this.P,1);
            this.P = CreateParamSet(this.P);
            try 
                this.P.pts = this.Sys.p;
                % Add property params
                props = this.Specs.values;
                for i=1:numel(props)
                    phi = props{i}
                    params_prop = get_params(phi);
                    this.SetParamSpec(fieldnames(params_prop)', cellfun(@(c) (params_prop.(c)), fieldnames(params_prop)),1);
                end
            end
            for ip = 1:numel(this.P.ParamList)
                this.Domains(ip).domain=[];  
            end
            this.resetStatus();
        end
        
        
        function ResetSimulations(this)
            % Removes computed trajectories
            this.P = SPurge(this.P);
            this.SignalRanges = [];
        end
        
        function ResetSelected(this)
            nb_pts = this.GetNbParamVectors();
            this.P.selected = zeros(1,nb_pts);
        end
        
        %%  Compare (FIXME)
        function cmp = compare(this, other)
            % Compares two BreachSet. Goes through a series of tests, logs
            % results and returns when an outstanding difference result is
            % found.
            % CHECKME: status number are quite arbitrary. Should evolve with
            % usage. Current guideline: large number means significant
            % (structural) difference, 0 means identical, large negative
            % means potiential bug, other can be anything, though structure
            % should be the same.
            
            cmp = BreachStatus;
            
            % checks if the handles are the same
            is_same_obj = (this==other);
            if is_same_obj
                cmp.addStatus(0, 'The two sets are the same objects in memory.')
                return;
            end
            
            % checks presence of P
            if (isempty(this.P))&&((isempty(other.P)))
                cmp.addStatus(0, 'The two objects have an empty parameter set.')
                return;
            elseif (isempty(this.P))
                cmp.addStatus(10000, 'Parameter set is empty in first object.');
                return;
            elseif (isempty(other.P))
                cmp.addStatus(10000, 'Parameter set is empty in second object.');
                return;
            end
            
            % Checks signals and parameters
            thisParamList = this.P.ParamList;
            otherParamList = other.P.ParamList;
            
            diff_pl = ~(isequal(thisParamList,otherParamList));
            if diff_pl
                % Checks signals
                sigthis = this.GetSignalNames();
                sigother = other.GetSignalNames();
                diff_signal = ~(isequal(sigthis,sigother));
                
                if diff_signal
                    cmp.addStatus(1, 'The two sets have different signals.')
                    return;
                else
                    thisParamList = this.P.ParamList(1:this.P.DimP);
                    otherParamList = other.P.ParamList(1:other.P.DimP);
                    diffParam = ~(isequal(thisParamList,otherParamList));
                    if diffParam
                        cmp.addStatus(1000, 'The two sets have different system parameter lists.')
                        return;
                    else
                        cmp.addStatus(10, 'The two sets have different property parameter lists.')
                    end
                end
                
            end
            
            % checks if parameter sets are equal
            if (isequal(this.P,other.P))
                cmp.addStatus(0,'The two sets parameters and traces are strictly identical.');
                return;
            end
            
            % checks pts nb
            nb_pts_this  = this.GetNbParamVectors();
            nb_pts_other = other.GetNbParamVectors();
            
            if nb_pts_this ~= nb_pts_other
                cmp.addStatus(10,'The two sets have different number of parameter vectors.');
                return;
            end
            
            % Checks system parameters
            rg_sys = this.P.DimX+1:this.P.DimP; % at this point, this is the same as other
            sys_pts_this = this.P.pts(rg_sys,:);
            sys_pts_other = other.P.pts(rg_sys,:);
            
            diff_sys_pts = norm(sys_pts_this-sys_pts_other);
            if (diff_sys_pts == 0)
                cmp.addStatus(0, 'The two sets have the same system parameter vectors.');
            else
                cmp.addStatus(1, ['Distance between the two system parameter vectors: ' num2str(diff_sys_pts)])
            end
            
            % Checks spec. parameters
            rg_pspec = this.P.DimP+1:size(this.P.pts,1); % at this point, this is the same as other
            if ~isempty(rg_pspec)
                spec_pts_this = this.P.pts(rg_pspec,:);
                spec_pts_other = other.P.pts(rg_pspec,:);
                
                diff_spec_pts = norm(spec_pts_this-spec_pts_other);
                if (diff_spec_pts == 0)
                    cmp.addStatus(0, 'The two objects have the same spec. parameter vectors.');
                else
                    cmp.addStatus(1, ['Distance between the two spec. parameter vectors: ' num2str(diff_spec_pts)])
                end
            end
            
            
            % Checks trajs field
            if (isfield(this.P, 'traj'))&&(~isfield(other.P, 'traj'))
                cmp.addStatus(1000, 'This has computed trajectories or traces but not other');
                return;
            elseif  (~isfield(this.P, 'traj'))&&(isfield(other.P, 'traj'))
                cmp.addStatus(1000, 'This has no computed trajectories or traces while other does have computed trajectories or traces.');
                return;
            elseif (isfield(this.P, 'traj'))&&(isfield(other.P, 'traj'))
                % both have traces, compare them
                nb_traj_this  = numel(this.P.traj);
                nb_traj_other = numel(other.P.traj);
                diff_nb_traj =nb_traj_this - nb_traj_other;
                if diff_nb_traj ~= 0
                    cmp.addStatus(100, 'Different numbers of traces.');
                    return;
                end
                
                if isequal(this.P.traj, other.P.traj)
                    cmp.addStatus(0, 'Traces are identical.');
                else
                    max_diff_param = 0;
                    max_diff_time = 0;
                    max_diff_X = 0;
                    for itraj=1:nb_traj_this
                        tr1 = this.P.traj{itraj};
                        tr2 = other.P.traj{itraj};
                        if (isequal(size(tr1.param), size(tr2.param)) &&  isequal(size(tr1.time), size(tr2.time)) && isequal(size(tr1.X), size(tr2.X)))
                            diff_p = tr1.param - tr2.param;
                            diff_time = tr1.time - tr2.time;
                            diff_X = tr1.X - tr2.X;
                            max_diff_param = max([max_diff_param, norm(diff_p)]);
                            max_diff_time = max([max_diff_time, norm(diff_time)]);
                            max_diff_X = max([max_diff_X, norm(diff_X)]);
                        else
                            cmp.addStatus(10, ['Traces ' itraj ' have different dimensions somehow.']);
                            return;
                        end
                    end
                    cmp.addStatus(1,['Max difference between trajectories: p:' num2str(max_diff_param) ' time:' num2str(max_diff_time) ' X:' num2str(max_diff_X)]);
                end
            end
            
            % Checks props fields
            
            if (isfield(this.P, 'props'))&&(~isfield(other.P, 'props'))
                cmp.addStatus(100, 'This has properties evaluated but not other');
                return;
            elseif  (~isfield(this.P, 'props'))&&(isfield(other.P, 'props'))
                cmp.addStatus(100, 'This has no properties evaluated while other does.');
                return;
            elseif (isfield(this.P, 'props'))&&(isfield(other.P, 'props'))
                % checks properties are the same
                if (isfield(this.P, 'props_names'))&&(isfield(other.P, 'props_names'))
                    if isequal(this.P.props_names,other.P.props_names)
                        if (isfield(this.P, 'props_values'))&&(isfield(other.P, 'props_values'))
                            
                            % both have traces, compare them
                            
                            this_props_vals  = this.P.props_values;
                            other_props_vals = other.P.props_values;
                            
                            if ~isequal(size(this_props_vals), size(other_props_vals))
                                cmp.addStatus(-100, 'Different numbers of property evaluations.');
                                return;
                            end
                            
                            if isequal(this_props_vals, other_props_vals)
                                cmp.addStatus(0, 'Property evaluations are identical.');
                            else
                                nb_traj = size(this_props_vals,2);
                                nb_phis = size(this_props_vals,1);
                                for itraj=1:nb_traj
                                    for iphi = 1:nb_phis
                                        tr1 = this_props_vals(iphi,itraj);
                                        tr2 = other_props_vals(iphi,itraj);
                                        if ~isequal(tr1,tr2)
                                            pre_status =  ['Trace: ' num2str(itraj) ' Property:' this.P.props_names{iphi}];
                                            if ~isequal(tr1.tau,tr2.tau)
                                                cmp.addStatus(1, [pre_status ' -- Evaluations at different times.']);
                                            elseif sign(tr1.val(1)) ~= sign(tr2.val(1))
                                                cmp.addStatus(1, [pre_status '-- Boolean satisfactions at origin are different']);
                                            elseif tr1.val(1) ~= tr2.val(1)
                                                cmp.addStatus(1, [pre_status '-- Quantitative satisfactions at origin are different']);
                                            end
                                            diff_val = norm(tr1.val-tr2.val);
                                            if diff_val
                                                cmp.addStatus(1, [pre_status '-- Overall difference in quantitative satisfactions: ' num2str(diff_val)]);
                                            end
                                            
                                        end
                                    end
                                end
                            end
                        else
                            cmp.addStatus(-100,'BUG: field props is present but not field props_values'); % probably dead code...
                            
                        end
                    else
                        cmp.addStatus(100,'Property names evaluated in this and other are different.');
                    end
                    
                else
                    cmp.addStatus(-100,'BUG: field props is present but not field props_names'); % probably dead code...
                end
            end
            
            
        end
        
    end
end
