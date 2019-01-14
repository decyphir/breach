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
    
    properties (SetObservable)
        P % legacy parameter structure - contains points data (in P.pts) and traces (P.traj) and many other fields whose purpose is slowly falling into oblivion
    end
    
    properties
        Domains = BreachDomain('double', [])
        ParamGens % optional classes performing parameter value transformation after a SetParam
        SignalRanges % ranges of values taken by each signal variable
        AppendWhenSample=false % when true, sampling appends new param vectors, otherwise replace.
        log_folder
        sigMap 
        sigMapInv
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
            
            this.sigMap = containers.Map();
            this.sigMapInv = containers.Map();
            
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
        function SetDomain(this, params, type, domain, enum)
            % SetDomain(params, breach_domains) or SetDomain(params, types, domains, enums)
            
            % First let's collect idx for params
            if nargin <3
                error('BreachSet:SetDomain:not_enough_args', 'SetDomain requires at least parameter names and one or several types.' );
            else
                
                params = check_arg_is_cell(params);
                idxs = zeros(1, numel(params));
                for ip = 1:numel(params)
                    param = params{ip};
                    [idx, found] = FindParam(this.P, param);
                    if found==0
                        error('BreachSet:SetDomain:param_or_signal_not_found', ['Parameter or signal '  param ' not found.']);
                    end
                    idxs(ip) = idx;
                end
            end
            
            % create domains
            switch nargin
                case 3
                    % is type a BreachDomain already?
                    if isa(type,'BreachDomain')||(iscell(type)&&~isempty(type)&&isa(type{1}, 'BreachDomain'))
                        if iscell(type)  % FIXME: this mess should be avoided by using cell everywhere for domains..
                            cell_type = type;   
                            type = BreachDomain;
                            for ic = 1:numel(cell_type) % converts cell of domains into array...
                                type(ic)= cell_type{ic};
                            end
                        end
                        
                        if numel(type)==1
                            type = repmat(type,1, numel(params));
                        end
                        if numel(type) ~= numel(params)
                            error('number of types and parameters or signals mismatch.')
                        end
                        
                        % Now we have as many BreachDomain objects as parameters
                        for ip = 1:numel(idxs)
                            this.Domains(idxs(ip)) = type(ip);
                        end
                    else  % type has to be a string or cell array
                        if ischar(type)
                            for ip = 1:numel(params)
                                this.Domains(idxs(ip)) = BreachDomain(type);
                            end
                        elseif isnumeric(type)
                            %  type given as ranges, go for double
                            if isequal(size(type), [1 2])
                                for ip = 1:numel(params)
                                    this.Domains(idxs(ip)) = BreachDomain('double', type);
                                end
                            elseif size(type, 2) == numel(params)
                                for ip = 1:numel(params)
                                    this.Domains(idxs(ip)) = BreachDomain('double', type(ip,:));
                                end
                            else
                                error('BreachSet:SetDomain:wrong_domain_size', 'Domain should be of size num_param x 2.');
                            end
                        end
                        
                        
                    end
                case 4
                    type = check_arg_is_cell(type, numel(params));
                    if isnumeric(domain)
                        %  domain given as ranges or enums
                        if isequal(size(domain), [1 2])
                            for ip = 1:numel(params)
                                this.Domains(idxs(ip)) = BreachDomain(type{ip}, domain);
                            end
                        elseif size(domain, 1) == numel(params)
                            for ip = 1:numel(params)
                                this.Domains(idxs(ip)) = BreachDomain(type{ip}, domain(ip,:));
                            end
                        else % might be enum, forward this to BreachDomain constructor
                            domain = check_arg_is_cell(domain, numel(params));
                            for ip = 1:numel(params)
                                this.Domains(idxs(ip)) = BreachDomain(type{ip}, domain{ip});
                            end
                        end
                    else
                        domain = check_arg_is_cell(domain, numel(params));
                        for ip = 1:numel(params)
                            this.Domains(idxs(ip)) = BreachDomain(type{ip}, domain{ip});
                        end
                    end
                    
                case 5
                    type = check_arg_is_cell(type, numel(params));
                    domain = check_arg_is_cell(domain, numel(params));
                    enum = check_arg_is_cell(enum, numel(params));
                    for ip = 1:numel(params)
                        this.Domains(idxs(ip)) = BreachDomain(type{ip}, domain{ip}, enum{ip});
                    end
            end
        end
        
        function dom = GetDomain(this, param)
            % BreachSet.GetDomain
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
            % BreachSet.CheckinDomain() Enforce parameters and signals to
            % adhere to their domains
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
            % BreachSet.CheckinDomainTraj()  Enforce signals to adhere to their domains
            
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
        
        %%  Params
        function SetParam(this, params, values, is_spec_param)
            % BreachSet.SetParam(params, values,  is_spec_param) sets values to
            % parameters listed in params. If the set contains only one sample,
            % creates as many sample as there are values. If the set has
            % several samples and there is only one value, set this value to
            % all samples. Otherwise, returns an error.
            
            if (~exist('is_spec_param', 'var'))
                is_spec_param = false;
            end
            
            ip = FindParam(this.P, params);
            i_not_sys = find(ip>this.P.DimP);
            if ~isempty(i_not_sys)
                if iscell(params)
                    nparam = params{i_not_sys(1)};
                else
                    nparam = params;
                end
                if isequal(is_spec_param,false)
                    warning('SetParam:param_not_in_list',['Parameter ' nparam ' was set but is not a system parameter.' ...
                        ' If this is intended to be a spec. parameter, consider using SetParamSpec instead.']);
                end
                for ii = 1:length(i_not_sys)
                    this.Domains(end+1) = BreachDomain(); % default domain for new parameters
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
            saved_traj = false;
            
            if ischar(is_spec_param)&&strcmp(is_spec_param, 'combine')
                if this.hasTraj()
                    P0 = this.P;
                    saved_traj = true;
                 end
                idx = N2Nn(2, [num_pts num_values]);
                old_pts = this.P.pts;
                this.P = Sselect(SPurge(this.P),1);
                this.P.pts = old_pts(:, idx(1,:));
                this.P.epsi= repmat(this.P.epsi,1, size(idx, 2));
                this.P.selected = zeros(1, size(idx, 2));
                this.P = SetParam(this.P, params, values(:, idx(2,:)));
            else  % legacy, i.e., not combine version
                if num_values==1 || num_values == num_pts
                    this.P = SetParam(this.P, params, values);
                elseif num_pts==1    % note in this case, we have to remove traces ( or see if maybe not, )
                    this.P = Sselect(SPurge(this.P),1);
                    this.P.pts = repmat(this.P.pts,1, size(values, 2));
                    this.P.epsi= repmat(this.P.epsi,1, size(values, 2));
                    this.P.selected = zeros(1, size(values, 2));
                    this.P = SetParam(this.P, params, values);
                else
                    error('SetParam:wrong_arguments_size', 'Dimension mismatch between values and parameters.');
                end
            end
            
            this.ApplyParamGens(params);
            
            % restore traj if needed
            if saved_traj
                this.P = Pfix_traj_ref(this.P, P0);
            end
            
        end
        
        function SetParamCfg(this, list_cfg)
            % SetParamCfg applies a cfg structure to set parameters. Struct
            % must have *char* fields params and values
            %
            if ~iscell(list_cfg)
                list_cfg = {list_cfg};
            end
            
            this.Reset();
            B0 = this.copy();
            % first operation
            cfg = list_cfg{1};
            B = B0.copy();
            for ip = 1:numel(cfg.params) % quick and dirty implementation, set parameters one by one using combine (grid) option
                p = cfg.params{ip};
                v = eval([ '[' cfg.values{ip} ']']);
                B.SetParam(p, v, 'combine');
            end
            
            % subsequent operations
            this.P= B.P;
            for ic = 2:numel(list_cfg)
                cfg = list_cfg{ic};
                B = B0.copy();
                for ip = 1:numel(cfg.params) % quick and dirty implementation, set parameters one by one using combine (grid) option
                    p = cfg.params{ip};
                    v = eval([ '[' cfg.values{ip} ']']);
                    B.SetParam(p, v, 'combine');
                end
                this.Concat(B);
            end
        end
        
        
        function SetParamSpec(this, params, values, ignore_sys_param)
            % BreachSet.SetParamSpec
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
           this.ApplyParamGens(params);
            
        end
 
        function SetParamGen(this, pg)
            if iscell(pg) 
                cellfun(@this.SetParamGenItem, pg);
            else 
                this.SetParamGenItem(pg);
            end
            
        end
        
        function SetParamGenItem(this, pg)
            this.ParamGens{end+1} = pg;
            
            % create/update domain of input parameters
            this.SetParam(pg.params, pg.p0, true);
            this.SetDomain(pg.params, pg.domain);
            
            % update domain of output parameters
            if ~isempty(pg.domain_out)
               for ip =1:numel(pg.params_out)
                   this.SetDomain(pg.params_out{ip}, pg.domain_out{ip});                   
               end
            end
                
            
            this.ApplyParamGens();
        end
        
        function ApplyParamGens(this, params)
            if ~isempty(this.ParamGens)
                if nargin==1
                    params = this.GetParamList(); 
                end
                
                % ensures params are names and ip are indices
                if ~isnumeric(params)
                    ip = FindParam(this.P, params);
                else
                    ip = params;
                    params = this.P.ParamList{ip};
                end
                for ig = 1:numel(this.ParamGens)
                    pg = this.ParamGens{ig};
                    params_in= pg.params;
                    if ~isempty(intersect(params, params_in))
                       p_in = this.GetParam(params_in);
                       p_out = pg.computeParams(p_in);
                       ip_out = FindParam(this.P, pg.params_out);
                       this.P.pts(ip_out,:) = p_out; 
                    end
                end
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
            % BreachSet.SetParamRanges set intervals for parameters (set domains as
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
            
            % kept for backward compatibility with legacy stuff
            this.P.dim = i_params;
            this.CheckinDomainParam();
        end
        
        function ranges = GetParamRanges(this, params)
            % BreachSet.GetParamRanges
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
        
        function params = GetParamList(this)
            % GetParamList returns parameter names
            params = this.P.ParamList(this.P.DimX+1:end);
        end
        
        function [params, idx] = GetSysParamList(this)
            % GetSysParamList 
            idx = this.P.DimX+1:this.P.DimP;
            params = this.P.ParamList(idx);
        end
        
        function prop_params = GetPropParamList(this)
            % GetSysParamList returns system parameter names
            prop_params = this.P.ParamList(this.P.DimP+1:end);
        end
        
        % Get the number of param vectors - -1 means P is empty
        function nb_pts = GetNbParamVectors(this)
            if isempty(this.P)
                nb_pts = -1;
            else
                nb_pts= size(this.P.pts,2);
            end
        end
        
        
        function [params, ipr] = GetVariables(this)
            [params, ipr] = GetBoundedDomains(this);
            if this.GetNbParamVectors()>1
                x = this.P.pts(this.P.DimX+1:end,:);
                ipr2 = find(sum(diff(x,1,2)~=0, 2)~=0)'+this.P.DimX;
                if ~isempty(ipr2)
                    ipr = union(ipr, ipr2, 'stable');
                    params = this.P.ParamList(ipr);
                end
            end
        end
        
        function [params, ipr] = GetSysVariables(this)
            [params, ipr] = GetVariables(this);
            if ~isempty(params)
                req_params = this.GetPropParamList();
                [params, i_diff] = setdiff(params, req_params);
                ipr = ipr(i_diff);
            end
        end
        
        function [params, ipr] = GetReqVariables(this)
            [params, ipr] = GetVariables(this);
            if ~isempty(params)
                req_params = this.GetPropParamList();
                [params, i_intersect] = intersect(params, req_params);
                ipr = ipr(i_intersect);
            end
        end
        
        function [ params, ipr]  = GetBoundedDomains(this)
            % GetNonEmptyDomains
            ipr = cellfun(@(c)(~isempty(c)), {this.Domains.domain});
            ipr = find(ipr);
            params =   this.P.ParamList(ipr);
        end
                
        %% Signals
        
        function  this =  SetSignalMap(this, varargin)
            % SetSignalMap defines aliases for signals - used by GetSignalValues 
            %
            % Input:  a map or a list of pairs or two cells
            
            arg_err_msg = 'Argument should be a containers.Map object, or a list of pairs of signal names, or two cells of signal names with same size.';
            switch nargin
                case 2
                    if ~isa(varargin{1}, 'containers.Map')
                        error('SetSignalMap:wrong_arg', arg_err_msg);
                    else
                        this.sigNamesMap = varargin{1};
                    end
                case 3
                    if iscell(varargin{2})
                        if ~iscell(varargin{2})||(numel(varargin{1}) ~= numel(varargin{2}))
                            error('SetSignalMap:wrong_arg', arg_err_msg);
                        end
                        for is = 1:numel(varargin{2})
                            if ~strcmp(varargin{1}{is},varargin{2}{is})
                                this.sigMap(varargin{1}{is}) = varargin{2}{is};
                                this.sigMapInv( varargin{2}{is} ) = varargin{1}{is};
                            end
                        end
                    else
                        if ischar(varargin{1})&&ischar(varargin{2})
                            if ~strcmp(varargin{1},varargin{2})
                                this.sigMap(varargin{1}) = varargin{2};
                                this.sigMapInv(varargin{2}) = varargin{1};
                            end
                        else
                            error('SetSignalMap:wrong_arg', arg_err_msg);
                        end
                    end
                otherwise
                    for is = 1:numel(varargin)/2
                        try
                            if ~strcmp(varargin{2*is-1},varargin{2*is})
                                this.sigMap(varargin{2*is-1}) = varargin{2*is};
                                this.sigMapInv(varargin{2*is}) = varargin{2*is-1};
                            end
                        catch
                            error('SetSignalMap:wrong_arg', arg_err_msg);
                        end
                    end
            end
     
            if this.verbose >= 2
                this.PrintSigMap();
            end
            
        end
        
        function ResetSigMap(this)
            this.sigMap = containers.Map();
        end

        function PrintSigMap(this)
            st = 'Signals Map:\n';
            
            keys  = this.sigMap.keys;
            for ip = 1:numel(keys)
                st =   sprintf([ st '%s ---> %s\n' ],keys{ip}, this.sigMap(keys{ip}));
            end
            fprintf(st);
            
        end
     
        function traces = GetTraces(this)
            % Get computed trajectories
            traces= [];
            if isfield(this.P,'traj')
                traces = this.P.traj;
            end
        end
        
        function [idx_ok, idx_sim_error, idx_invalid_input, st_status]  = GetTraceStatus(this)
            % BreachSet.GetTraceStatus returns indices of ok traces, error and input invalid.
            idx_ok = [];
            idx_sim_error = [];
            idx_invalid_input = [];
            
            if this.hasTraj
                nb_pts = size(this.P.pts, 2);
                idx_ok = 1:nb_pts;
                if isfield(this.P.traj{1}, 'status')
                    
                    traj_status = zeros(1, nb_pts);
                    for it = 1:nb_pts
                        if this.P.traj_ref(it)
                            traj_status(it) = this.P.traj{this.P.traj_ref(it)}.status;
                        end
                    end
                    idx_ok = find(traj_status == 0);
                    idx_sim_error = find(traj_status == -1);
                    idx_invalid_input = find(traj_status == -2);
                end
            end
            
            st_status = '';
            if ~isempty(idx_sim_error)
                st_status = sprintf('%s %d samples caused a simulation error.', st_status, numel(idx_sim_error));
            end
            if ~isempty(idx_invalid_input)
                st_status = sprintf('%s %d simulations skipped for invalid inputs.', st_status, numel(idx_invalid_input));
            end
            
        end
        
        function dispTraceStatus(this)
            [~,~, ~, st_status]  = GetTraceStatus(this);
            if ~isempty(st_status)
                fprintf('%s\n', st_status);
            end
        end
        
        function [Bok, Bsim_error, Binvalid_input] = FilterTraceStatus(this)
            % BreachSet.FilterTraceStatus() Removes and extract traces with simulator errors or input not satisfying constraints
            %
            
            Bok = [];
            Bsim_error =[];
            Binvalid_input = [];
            [idx_ok, idx_sim_error, idx_invalid_input]  = GetTraceStatus(this);
            
            if ~isempty(idx_ok)
                Bok = this.ExtractSubset(idx_ok);
            end
            
            if ~isempty(idx_sim_error)
                Bsim_error = this.ExtractSubset(idx_sim_error);
            end
            
            if ~isempty(idx_invalid_input)
                Binvalid_input = this.ExtractSubset(idx_invalid_input);
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
            SigNames = this.GetSignalList();
        end
        
        function [signals, idx] = GetSignalList(this)
            % GetSignalList returns signal names
            signals = this.P.ParamList(1:this.P.DimX);
            idx = 1:this.P.DimX;
        end
        
        function SigNames = GetAllSignalsList(this)
        % GetAllSignalsList returns all signals names including aliases    
            SigNames = this.expand_signal_name('.*');
        end
            
        function X = GetSignalValues(this, signals, itrajs, t)
            % BreachSet.GetSignalValues(signals, idx_traces, time) - in case of several trajectories, return cell array
            if (~isfield(this.P,'traj'))
                error('GetTrajValues:NoTrajField','Compute/import trajectories first.')
            end
            
            if ischar(signals)
                signals = {signals};
            end
            
            if iscell(signals)
                [signals_idx, type] = this.FindSignalsIdx(signals);
                if any(type==0)
                    not_found= find(type==0);
                    if ischar(signals)
                        sig = signals;
                    else
                        sig = signals{not_found(1)};
                    end
                    error('GetSignalValues:sig_not_found', 'Signal %s not found.',sig);
                end
            elseif isnumeric(signals)
                signals_idx=signals;
            else
                error('GetSignalValues:signals_type', 'signals should be a string or cell of string or numeric array ')
            end
            
            if ~exist('itrajs','var')
                itrajs= 1:numel(this.P.traj);
            end
            
            nb_traj = numel(itrajs);
            
            X = cell(nb_traj,1);
            for i_traj = 1:numel(itrajs)
                Xi = this.P.traj{itrajs(i_traj)}.X;
                if (~exist('t','var'))
                    X{i_traj} = Xi(signals_idx,:);
                else
                    X{i_traj} = interp1(this.P.traj{itrajs(i_traj)}.time, Xi(signals_idx,:)',t)';
                    if numel(signals_idx)==1
                        X{i_traj} = X{i_traj}';
                    end
                end
            end
            if nb_traj==1
                X = X{1};
            end
        end
        
        function [idx, ifound] = FindSignalsIdx(this, signals)
            % resolve sigMap
            if ischar(signals)
                signals = {signals};
            end
            % For all signals, first use sigMap-less search, then check
            % aliases
            for isig = 1:numel(signals)
                sig = signals{isig}; 
                [idx(isig), ifound(isig)] = FindParam(this.P, sig);
                aliases_sig = setdiff(this.getAliases(sig), sig);
                for ais = 1:numel(aliases_sig)
                    [idx_s, ifound_s] = FindParam(this.P, aliases_sig{ais});
                    if ifound_s
                        ifound(isig)=true;
                        idx(isig)=idx_s;
                    end
                end
            end
        end
        
        function SaveSignals(this, signals, folder, name,i_trajs)
            % BreachSet SaveSignals Save signals in mat files as simple time series
            
            if ~this.hasTraj()
                error('No signals computed for this set.' );
            end
            
            % arguments
            if ~exist('signals','var')||isempty(signals)
                signals = this.GetSignalList();
            end
            if ~exist('folder','var')||isempty(folder)
                folder =  [this.whoamI '_Signals_' datestr(now, 'dd_mm_yyyy_HHMM')];
            end
            if ~exist('name','var')||isempty(name)
                name= 'Signals_';
            end
            if ~exist('i_trajs','var')||isempty(i_trajs)
                i_trajs = 1:numel(this.P.traj);
            end
            
            [success,msg,msg_id] = mkdir(folder);
            if success == 1
                if isequal(msg_id, 'MATLAB:MKDIR:DirectoryExists')
                    this.disp_msg(['Saving in existing folder ' folder]);
                else
                    this.disp_msg(['Saving in new folder ' folder]);
                end
            else
                error(['Couldn''t create folder'  folder '. mkdir returned error: ' msg]);
            end
            
            [sys_param_list, sys_param_idx] = this.GetSysParamList();
            for ip = 1:numel(i_trajs)
                fname = [folder filesep name num2str(ip) '.mat'];
                time = this.P.traj{ip}.time';
                save(fname, 'time');
                for is = 1:numel(signals)
                    sig = this.GetSignalValues(signals{is}, ip)';
                    eval([ signals{is} '= sig;']);
                    save(fname,'-append',  signals{is});
                end
                % save parameter values
                for iparam=1:numel(sys_param_idx)
                    p = this.P.traj{ip}.param(sys_param_idx(iparam));
                    if ~isequal(sys_param_list{iparam}, 'file_idx')
                        eval([ sys_param_list{iparam} '= p;']);
                        save(fname,'-append',  sys_param_list{iparam});
                    end
                end
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
        function SampleDomain(this, params, num_samples, method, opt_multi, max_num_samples)
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
            elseif isequal(method, 'corners')
                if exist('max_num_samples', 'var')
                    x = sample(domains{:}, num_samples, method, max_num_samples);
                else
                    x = sample(domains{:}, num_samples, method);
                end
                
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
                    num_new = size(x,2);
                    if exist('max_num_samples', 'var')
                        idx = N2Nn(2, [num_old num_new],max_num_samples);
                    else
                        idx = N2Nn(2, [num_old num_new]);
                    end
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
        function CornerSample(this, max_num_samples)
            
            if nargin ==1
                max_num_samples = inf;
            end
            bnd_params = this.GetBoundedDomains();
            if this.AppendWhenSample
                this.SampleDomain(bnd_params, 2, 'corners', 'append', max_num_samples);
            else
                this.SampleDomain(bnd_params,2, 'corners', 'replace', max_num_samples);
            end
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
                
        %% Concatenation, ExtractSubset - needs some additional compatibility checks...
        function Concat(this, other)
            this.P = SConcat(this.P, other.P);
        end
        
        function other  = ExtractSubset(this, idx)
            other = this.copy();
            other.P = Sselect(this.P, idx);
        end
        
        %% Plot parameters
        function PlotParams(this, varargin)
            % Plot parameters
            gca;
            %P = DiscrimPropValues(this.P);
            params = SplotPts(this.P, varargin{:});
            
            %% Datacursor mode customization
            h = datacursormode(gcf);
            h.UpdateFcn = @myupdatefcn;
            h.SnapToDataVertex = 'on';
            datacursormode on
            
            function [txt] = myupdatefcn(obj,event_obj)
                
                pos = event_obj.Position;
                switch numel(params)
                    case 1
                        txt = {[params{1} ':' num2str(pos(1))],...
                            };
                    case 2
                        txt = {[params{1} ':' num2str(pos(1))],...
                            [params{2} ': ',num2str(pos(2))],...
                            };
                        
                    case 3
                        txt = {[params{1} ':' num2str(pos(1))],...
                            [params{2} ': ',num2str(pos(2))],...
                            [params{3} ': ',num2str(pos(3))],...
                            };
                end
            end
            
        end
        
        %% Plot domains
        function PlotDomain(this, params)
            % BreachSet.PlotMixedDomain UNFINISHED, tries to plot enum/int
            % domain differently from dense domains
            
            gca;
            set(gca, 'XLimMode', 'auto', 'YLimMode', 'auto', 'ZLimMode', 'auto')
            
            % default style
            col = [0 0 1];
            alpha = 0.05;
            
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
                    
                    %% unzoom
                    xlim = get(gca, 'XLim');
                    dxlim = xlim(2) - xlim(1);
                    set(gca, 'XLim', [xlim(1)-dxlim/30, xlim(2)+dxlim/30]);
                case 2 % two domains
                    domain1 = this.GetDomain(params{1});
                    domain2 = this.GetDomain(params{2});
                    
                    start = [domain1.domain(1), domain2.domain(1)];
                    sz = [domain1.domain(2) - domain1.domain(1),domain2.domain(2) - domain2.domain(1)];
                    d= rect(start, sz, col,alpha);
                    
                    xlabel(params{1}, 'Interpreter', 'none');
                    ylabel(params{2}, 'Interpreter', 'none');
                    grid on;
                    
                    %% Unzoom slightly
                    xlim = get(gca, 'XLim');
                    dxlim = xlim(2) - xlim(1);
                    set(gca, 'XLim', [xlim(1)-dxlim/30, xlim(2)+dxlim/30]);
                    ylim = get(gca, 'YLim');
                    dylim = ylim(2) - ylim(1);
                    set(gca, 'YLim', [ylim(1)-dylim/30, ylim(2)+dylim/30]);
                    
                otherwise
                    
                    domain1 = this.GetDomain(params{1});
                    domain2 = this.GetDomain(params{2});
                    domain3 = this.GetDomain(params{3});
                    
                    start = [domain1.domain(1), domain2.domain(1), domain3.domain(1)];
                    sz = [domain1.domain(2) - domain1.domain(1),domain2.domain(2) - domain2.domain(1),domain3.domain(2) - domain3.domain(1)];
                    d= voxel(start, sz, col, alpha);
                    
                    xlabel(params{1}, 'Interpreter', 'none');
                    ylabel(params{2}, 'Interpreter', 'none');
                    zlabel(params{3}, 'Interpreter', 'none');
                    grid on;
                    
                    xlim = get(gca, 'XLim');
                    dxlim = xlim(2) - xlim(1);
                    set(gca, 'XLim', [xlim(1)-dxlim/30, xlim(2)+dxlim/30]);
                    ylim = get(gca, 'YLim');
                    dylim = ylim(2) - ylim(1);
                    set(gca, 'YLim', [ylim(1)-dylim/30, ylim(2)+dylim/30]);
                    zlim = get(gca, 'ZLim');
                    dzlim = zlim(2) - zlim(1);
                    set(gca, 'ZLim', [zlim(1)-dzlim/30, zlim(2)+dzlim/30]);
                    view(-37.5, 30);
            end
            set(d,'EdgeAlpha',0.1);
            
            
            % 1d x direction
            function  plotxdomain(dom, y0)
                width_y = .1;
                start = [dom.domain(1), y0-width_y/2];
                sz = [dom.domain(2) - dom.domain(1), width_y];
                d= rect(start, sz, col, alpha);
                set(gca, 'YLim', [y0-10*width_y, y0+10*width_y]);
            end
            
        end
        
        function PlotMixedDomain(this, params)
            % BreachSet.PlotMixedDomain UNFINISHED, tries to plot enum/int
            % domain differently from dense domains
            
            gca;
            
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
                    view(-37.5, 30);
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
        
        %% Requirements       
        function  SortbyRob(this)
            sat_values = this.GetSatValues();
            [ ~, order_rob] = sort(sum(sat_values,1));
            this.P = Sselect(this.P, order_rob);
        end
        
        function  SortbySat(this)
            sat_values = this.GetSatValues();
            [ ~, order_rob] = sort(sum(sat_values>=0,1));
            this.P = Sselect(this.P, order_rob);
        end
        
        %% Printing/Exporting
        function [success, msg, msg_id] = SaveResults(this, folder_name, varargin)
            % FIXME does not support attributes? 
            
            % Additional options
            if ~exist('folder_name', 'var')
                folder_name = '';
            end
            options = struct('FolderName', folder_name, 'SaveBreachSystem', true, 'ExportToExcel', false, 'ExcelFileName', 'Results.xlsx');
            options = varargin2struct(options, varargin{:});
            
            if isempty(options.FolderName)
                try
                    options.FolderName = [this.mdl.name '_Results_' datestr(now, 'dd_mm_yyyy_HHMM')];
                catch
                    options.FolderName = [this.whoamI '_Results_' datestr(now, 'dd_mm_yyyy_HHMM')];
                end
            end
            
            folder_name = options.FolderName;
            [success,msg,msg_id] = mkdir(folder_name);
            trace_folder_name = [folder_name filesep 'traces'];
            [success,msg,msg_id] = mkdir(trace_folder_name);
            
            if success == 1
                if isequal(msg_id, 'MATLAB:MKDIR:DirectoryExists')
                    this.disp_msg(['Saving in existing result folder at ' folder_name]);
                else
                    this.disp_msg(['Created result folder at ' folder_name]);
                end
                this.log_folder = folder_name;
            else
                error(['Couldn''t create folder'  folder_name '.']);
            end
            
            if ~this.hasTraj()
                error('Breach:SaveResult:no_trace', 'No trace to save - run Sim command first');
                return;
            end
            
            [summary, traces] = this.ExportTracesToStruct();
            %saving summary
            summary_filename = [folder_name filesep 'summary'];
            save(summary_filename,'-struct', 'summary');
            
            if  options.SaveBreachSystem
                breachsys_filename  = [folder_name filesep 'breach_system'];
                breachsys_name = this.whoamI;
                evalin('base', ['save(''' breachsys_filename ''', ''' breachsys_name  ''', ''-v7.3'');'] );
            end
            
            for it=1:numel(traces)
                trace_filename = [trace_folder_name filesep num2str(it) '.mat'];
                trace = traces(it);
                save(trace_filename,'-struct', 'trace');
            end
            
            if options.ExportToExcel
                this.ExportToExcel(options.ExcelFileName);
            end
        end
                
        function [summary, traces] = ExportTracesToStruct(this,i_traces, varargin)
            % BreachSet.ExportTracesToStruct
            
            summary = this.GetSummary();
            traces = [];
            if ~this.hasTraj()
                error('Breach:ExportTrace:no_trace', 'No trace to export - run Sim command first');
                return;
            end
            
            num_traces = numel(this.P.traj);
            if nargin==1
                i_traces = 1:num_traces;
            end
            
            % Additional options
            options = struct('FolderName', '','IncludesOnlySignals', [], 'ExcludeSignals', []);
            options = varargin2struct(options, varargin{:});
            
            if isempty(options.FolderName)
                options.FolderName = ['Results_' datestr(now, 'dd_mm_yyyy_HHMM')];
            end
            
            %% Common stuff
            
            % parameter names
            param_names = this.GetSysParamList();
            
            % input signal names
            signal_names = this.GetSignalList();
            
            if isfield(this.P,'props_names')
                spec_names = this.P.props_names;
            end
            
            
            %% traces
            summary.filenames = {};
            summary.paths = {};
            for it = i_traces
                
                % check status
                if isfield(this.P.traj{it}, 'status')&&(this.P.traj{it}.status~=0)
                    warning('SaveResults:suspicious_trace','Trace %d has suspicious status, likely resulting from simulation error.', it)  
                end
                
                % params
                traces(it).params.names = param_names;
                traces(it).params.values = this.GetParam(param_names,it)';
                
                % time
                traces(it).time = this.P.traj{it}.time;
                
                % input signals
                traces(it).signals.names = signal_names;
                traces(it).signals.values =  this.GetSignalValues(signal_names, it);
                
                % specifications
                if isfield(this.P,'props_names')
                    traces(it).specs.ids = spec_names;
                    for ip = 1:numel(this.P.props_names)
                        traces(it).specs.stl_formula{ip} = disp(this.P.props(ip));
                        traces(it).specs.stl_formula_full{ip} = disp(this.P.props(ip),0);
                        params = get_params(this.P.props(ip));
                        traces(it).specs.params(ip).names = fieldnames(params);
                        traces(it).specs.params(ip).values = this.GetParam(fieldnames(params), it)';
                        traces(it).specs.rob(ip).time =this.P.props_values(ip, it).tau;
                        traces(it).specs.rob(ip).values =  this.P.props_values(ip, it).val;
                        traces(it).specs.status(ip) =  this.P.props_values(ip, it).val(1)>=0;
                    end
                end
            end
            
        end

        function traces = ExportTraces(this, signals, params, varargin)
            
            if ~exist('signals','var')
                signals = {}; % means all
            end
            if ~exist('params','var')||isempty(params)
                params = {}; % means all
            end
            
            % Options
            options = struct('WriteToFolder','');
            options = varargin2struct(options, varargin{:});
            
            if ~isempty(options.WriteToFolder)
                if ~exist(options.WriteToFolder,'dir' )
                    [success, err_msg] = mkdir(options.WriteToFolder);
                    if ~success
                        error('Folder creation failed with error:\n %s', err_msg);
                    end
                end
                dir_traces = options.WriteToFolder;
            else
                dir_traces = '';
            end
            
            [signature,~, params] = this.GetSignature(signals, params);
            num_traces = numel(this.P.traj);
            signals = signature.signals_reps; % signal representants, assuming there are aliases
            param_values = this.GetParam(params);
            for it = 1:num_traces
                traj = this.P.traj{it};
                X = this.GetSignalValues(signals, it);
                
                if ~isempty(dir_traces)
                    traces{it} = matfile([dir_traces filesep num2str(it) '.mat'], 'Writable',true);
                end
                
                if isfield(traj, 'status')
                    traces{it}.status = traj.status;
                end
                traces{it}.signature = signature;
                traces{it}.param = [zeros(1,numel(signals)) param_values(:,it)'];
                traces{it}.time = traj.time;
                traces{it}.X = X;
                
                if ~isempty(dir_traces)
                    traces{it}.Properties.Writable= false;
                end
            end
        end
        
        function [signature, signals, params] = GetSignature(this, signal_list, param_list)
            %  GetSignature returns information about signals and parameters 
            
            % gets signals signature
            if ~exist('signal_list', 'var')||isempty(signal_list)
                signal_list = this.GetSignalList();
            end
            
            [signals, unknown] = this.expand_signal_name(signal_list);
            if ~isempty(unknown)
                warning('GetSignature:signal_unknown', 'Signal or attribute %s not found.', unknown{1});
            end
            signature = this.GetSignalSignature(signals);
            
            % gets params signature
            if ~exist('param_list', 'var')||isempty(param_list)
                param_list = this.GetParamList();
            end

            [params, unknown] = this.expand_param_name(param_list);
            if ~isempty(unknown)
                warning('GetSignature:param_unknown', 'Parameter %s unknown.', unknown{1});
            end
            sigp = this.GetParamSignature(params);
            
            f = fieldnames(sigp);
            for i = 1:length(f)
                signature.(f{i}) = sigp.(f{i});
            end
        end
        
        function [sigs] = GetSignalSignature(this, signals)
           if nargin <=1
                signals = this.GetSignalList(); 
           end
                
           if ischar(signals)
                signals= {signals};
            end
            sigs.signals = signals;
            sigs.signal_attributes = {};
            sigs.signals_map_idx = [];
            sigs.signals_reps={};  
            for is = 1:numel(signals)
                sig= signals{is};
                sig_aliases = this.getAliases(sig);
                idx = inf; 
                for ia = 1:numel(sig_aliases)
                    idx_ia =  find(strcmp(sig_aliases(ia), signals),1);
                    if ~isempty(idx_ia)
                        idx = min(idx,idx_ia ); % find first position in signals an alias appears
                    end
                end
                if idx==is % first time we see this guy, take it as rep
                    sigs.signals_reps{end+1} = sig;
                    idx = numel(sigs.signals_reps);
                end
                
                if idx==inf
                    warning('BreachSet:GetSignalSignature:not_found', 'Signal or alias %s not found.', sig);
                end
                
                sigs.signals_map_idx(is) = idx; 
                dom = this.GetDomain(signals{idx});
                sigs.signal_types{is}  = dom.type;
                if isempty(dom.type)  % ? 
                    sigs.signal_types{is}  = 'double'; 
                end
                
                %  Add attributes indexes
                atts = this.get_signal_attributes(sig);
                sigs.signal_attributes = union(sigs.signal_attributes, atts);
                for ia = 1:numel(atts)
                    f = [atts{ia} 's_idx'];
                    if ~isfield(sigs, f)
                        sigs.(f) = is;
                    else
                        sigs.(f)(end+1) = is;
                    end
                end
            end
            if size(sigs.signal_attributes,1)>1
                sigs.signal_attributes = sigs.signal_attributes';
            end
        end
        
        function sigp = GetParamSignature(this, params)
            
            if nargin<=1
                params = this.GetParamList();
            end
            sigp.params=params;
            sigp.param_attributes ={};
            for ip = 1:numel(params)
                par  =params(ip);
                dom = this.GetDomain(par);
                sigp.param_types{ip}  = dom.type ;
                if isempty(dom.type)  % ? 
                    sigp.param_types{is}  = 'double'; 
                end
                %  Add attributes indexes
                atts = this.get_param_attributes(par);
                sigp.param_attributes = union(sigp.param_attributes, atts); 
                for ia = 1:numel(atts)
                    f = [atts{ia} 's_idx'];
                    if ~isfield(sigp, f)
                        sigp.(f) = ip;
                    else
                        sigp.(f)(end+1) = ip;
                    end
                end
                
            end
        end
        
        
        function summary = GetSummary(this)
            
            if this.hasTraj()
                num_traces = numel(this.P.traj);
            else
                num_traces = 0;
            end
            
            % parameter names
            param_names = this.GetSysParamList();
            
            % input signal names
            signal_names = this.GetSignalList();
            
            if isfield(this.P,'props_names')
                spec_names = this.P.props_names;
            end
            summary.signature= this.GetSignature();
            summary.date = datestr(now);
            summary.num_traces = num_traces;
            summary.params.names = param_names;
            summary.params.values = this.GetParam(summary.params.names);
        
            if isfield(this.P, 'props')
                summary.specs.names = spec_names;
                this.SortbyRob();
                this.SortbySat();
                summary.specs.rob = this.GetSatValues();
                summary.specs.sat = summary.specs.rob>=0;
                summary.num_sat = - sum( ~summary.specs.sat, 1  );
            end
        
        end
        
        function b = is_variable(this, par)
            b = false ;
            [i, f] = FindParam(this.P, par);
            if f
                dom = this.Domains(i);
                b = ~isempty(dom.domain);
                if ~b
                    b = any(diff(this.P.pts(i,:)));
                end
            end
        end
        
        function att = get_signal_attributes(this, sig)
            % returns nature to be included in signature
            att = {};
        end

        function atts = get_param_attributes(this, param)
            % returns nature to be included in signature
            atts = {};
            if this.is_variable(param)
                atts = [atts {'variable'}];
            else
               atts = [atts {'const'}];
            end
        end

        function [sig_names, unknown] =  expand_signal_name(this, signals)
            % expand_signal_name expands a string into a set of signal names by attribute or regular expression search
            sig_names = {};
            unknown = {};
            if ischar(signals)
                signals = {signals};
            end
            S = this.GetSignalSignature();
            
            for isig = 1:numel(signals)
                sig = signals{isig};
                if ismember(sig,S.signals)  % is a signal name already
                    sig_names = [sig_names {sig}];
                elseif ismember(sig, S.signal_attributes)  % attribute
                    sig_names = [sig_names S.signals(S.([sig 's_idx'])) ];
                else  % regexp search
                    sig_names = [sig_names S.signals(cellfun(@(c)(~isempty(c)),  regexp(S.signals,sig)))];
                    if isempty(sig_names)
                        unknown = [unknown sig];
                    end
                end
            end
            sig_names = unique(sig_names, 'stable');
            % adds in aliases
            sig_names = union(sig_names,this.getAliases(sig_names), 'stable');
            
        end
        
        function [param_names, unknown] =  expand_param_name(this, params)
            % expand_signal_name expands a string into a set of signal names by attribute or regular expression search
            param_names = {};
            unknown = {};
            if ischar(params)
                params = {params};
            end
            S = this.GetParamSignature();
            
            for ipar = 1:numel(params)
                param = params{ipar};
                if ismember(param,S.params)  % is a signal name already
                    param_names = [param_names {param}];
                elseif ismember(param, S.param_attributes)  % attribute
                    param_names = [param_names S.params(S.([param 's_idx'])) ];
                else  % regexp search
                    param_names = [param_names S.params(cellfun(@(c)(~isempty(c)),  regexp(S.params,param)))];
                    if isempty(param_names)
                        unknown = [unknown param];
                    end
                end
            end
        end
        
        function ExportToExcel(this, excel_file)
            [summary, traces] = this.ExportTracesToStruct();
            global BreachGlobOpt
            
            if ~isfield(summary, 'specs')
                warning('Breach:SaveResult:no_spec_for_Excel','Export to Excel requested but there is no requirement result to report. Excel file not created.');
            else
                breach_dir = BreachGlobOpt.breach_dir;
                template_file_path = [breach_dir filesep 'Ext' filesep 'Toolboxes' filesep 'ExportResults' filesep 'BreachResults_template.xlsx'];
                copyfile(template_file_path, excel_file);
                
                % Write header
                for ispec = 1:numel(summary.specs.names)
                    hdr{ispec} = ['Req. ' num2str(ispec)];
                end
                
                for iparam = ispec+1:ispec+numel(summary.test_params.names)
                    hdr{iparam} = ['param. ' num2str(iparam-ispec) ];
                end
                xlswrite(excel_file, hdr, 1, 'B1');
                xlswrite(excel_file, [summary.specs.names summary.test_params.names], 1, 'B2');
                
                % Write data
                xlswrite(excel_file, [ summary.num_sat' summary.specs.rob' summary.test_params.values'] , 1, 'A3');
                this.disp_msg(['Summary written into ' excel_file]);
            end
        end
        
        function st = PrintSignals(this)
            st = sprintf('---- SIGNALS ----\n');
            for isig = 1:this.P.DimX
                st = sprintf([st '%s %s\n'], this.P.ParamList{isig}, this.get_signal_attributes_string(this.P.ParamList{isig}));
            end
            st = sprintf([st '\n']);
            st = [st this.PrintAliases()];
            if nargout==0
                fprintf(st);
            end
            
        end
        
        function st = PrintAliases(this)
            st = '';
            if ~isempty(this.sigMap)
                st = sprintf('---- ALIASES ----\n');
                keys = union(this.sigMap.keys(), this.sigMapInv.keys());
                printed = {};
                for ik = 1:numel(keys)
                    if this.sigMap.isKey(keys{ik})
                        sig =  this.sigMap(keys{ik});
                    else
                        sig =this.sigMapInv(keys{ik});
                    end
                    
                     [idx, found] = this.FindSignalsIdx(sig);
                     if found
                         sig = this.P.ParamList{idx};
                         aliases = setdiff(this.getAliases(sig),sig);
                         al_st = cell2mat(cellfun(@(c) ([ c ', ']), aliases, 'UniformOutput', false));
                         al_st = al_st(1:end-2);
                     else
                         aliases = setdiff(this.getAliases(sig),sig);
                         al_st = cell2mat(cellfun(@(c) ([ c ', ']), aliases, 'UniformOutput', false));
                         al_st = al_st(1:end-2);
                         al_st = [al_st(1:end-2) ' (not linked to data)' ];
                     end
                     if ~ismember(sig, printed)
                        st = sprintf([st '%s <--> %s\n'], sig, al_st );
                        printed = [printed {sig} aliases];
                    end
                 end
                st = sprintf([st '\n']);
                if nargout==0
                    fprintf(st);
                end
            
            end
        end
        
        
        function st = get_signal_attributes_string(this, sig) 
            atts = this.get_signal_attributes(sig);
            if isempty(atts)
                st = '';
            else
                st  = '(';
                for ia = 1:numel(atts)-1
                    st = sprintf([st '%s,'], atts{ia});
                end
                st = sprintf([st '%s)'], atts{end});
            end
        end
        
        function st = PrintParams(this)
            st = '';
            nb_pts= this.GetNbParamVectors();
            if (nb_pts<=1)
                st = sprintf('-- PARAMETERS --\n');
                for ip = this.P.DimX+1:numel(this.P.ParamList)
                    st = sprintf([st '%s=%g       %s\n'],this.P.ParamList{ip},this.P.pts(ip,1), this.Domains(ip).short_disp(1));
                end
            else
                st = sprintf([st '-- PARAMETERS -- (%d vectors):\n'],nb_pts);
                for ip = this.P.DimX+1:numel(this.P.ParamList)
                    st = sprintf([st '%s     %s\n'],this.P.ParamList{ip}, this.Domains(ip).short_disp(1));
                end
            end
            
            st = sprintf([st ' \n']);
            
            if nargout==0
                fprintf(st);
            end
        end
        
        %% Misc
        function s= isSignal(this,params)
            idx_s = FindParam(this.P, params);
            s = idx_s <= this.P.DimX;
        end
        
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
        
        %% Reset functions - will need some cleaning and rationalizing eventually
        function Reset(this)
            % BreachSet.Reset()
            
            this.P = Sselect(this.P,1);
            this.P = CreateParamSet(this.P);
            try
                this.P.pts = this.Sys.p;
                % Add property params
                props = this.Specs.values;
                for i=1:numel(props)
                    phi = props{i};
                    params_prop = get_params(phi);
                    this.SetParamSpec(fieldnames(params_prop)', cellfun(@(c) (params_prop.(c)), fieldnames(params_prop)),1);
                end
            end
            for ip = 1:numel(this.P.ParamList)
                this.Domains(ip).domain=[];
            end
            
            this.resetStatus();
        end
        
        function ResetEpsi(this)
            % (Legacy) Set param ranges around individual parameter vectors to zero
            this.P.epsi(:,:) = 0;
        end
        
        function ResetDomains(this)
            % BreachSet.ResetDomains Sets all domains to empty
            for id = 1:numel(this.Domains)
                this.Domains(id).domain = [];
            end
        end
        
        function ResetDomain(this, params)
            % BreachSet.ResetDomains Sets given domain(s) to empty
            
            idx_params = FindParam(this.P, params);
            for id = 1:numel(idx_params)
                if (idx_params(id) <= numel(this.Domains)) % param was found, otherwise do nothing
                    this.Domains(idx_params(id)).domain = [];
                end
            end
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
        
        function aliases = getAliases(this, signals)
            if ischar(signals)
                signals = {signals};
            end
            
            aliases = signals;
            sig_queue = signals;
            
            while ~isempty(sig_queue)
                sig = sig_queue{1};
                sig_queue = sig_queue(2:end);
                if this.sigMap.isKey(sig)
                    nu_sig = this.sigMap(sig);
                    check_nusig()
                end
                if this.sigMapInv.isKey(sig)
                    nu_sig = this.sigMapInv(sig);
                    check_nusig()
                end
            end
            
            % check sigMapInv for double alias
            for  invkey = this.sigMapInv.keys()
                sig = this.sigMapInv(invkey{1});
                if ismember(sig,aliases)
                    aliases = union(aliases, invkey{1});
                end
            end
            
            function check_nusig()
                if ~ismember(nu_sig, aliases)
                    aliases = [aliases {nu_sig}];
                    sig_queue = [sig_queue nu_sig];
                end
            end
        end
         
    end
    methods (Access=protected)    
        
        
        
        function Xp = get_signals_from_traj(this, traj, names)
            idx = FindParam(this.P, names); %  not fool proof, but not supposed to be used by fools
            Xp = traj.X(idx,:);
        end
  
        function traj = set_signals_in_traj(this, traj, names, Xp)
            idx = FindParam(this.P, names); %  not fool proof, but not supposed to be used by fools
            traj.X(idx,:) = Xp;
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
