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
    %   AppendWhenSample=false   - when true, sampling appends new param samples, otherwise replace.
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
        AliasMap
    end

    properties  %  coverage stuff
        CoverOpts   %  options for coverage computation
        CoverRes     %   results for previous coverage computation
        DeltaGridMapObj
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
            this.AliasMap = containers.Map();

            this.DeltaGridMapObj = containers.Map();

            if nargin>=1 && ischar(Sys)
                Sys= {Sys};
            end


            switch nargin
                case 0
                    return;
                case 1
                    if isaSys(Sys)
                        this.P = CreateParamSet(Sys);
                    elseif iscell(Sys) % assumes parameter names
                        Sys = CreateSystem({}, Sys, zeros(1, numel(Sys)));
                        this.P = CreateParamSet(Sys);
                    end
                case 2
                    if isaSys(Sys)
                        this.P = CreateParamSet(Sys,params);
                    elseif iscell(Sys) % assumes parameter names
                        ranges = params;
                        params = Sys;
                        Sys = CreateSystem({}, Sys, zeros(1, numel(Sys)));
                        this.P = CreateParamSet(Sys,params,ranges);
                        this.SetParamRanges(params,ranges);
                    end
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
                        warning('BreachSet:SetDomain:param_or_signal_not_found', ['Parameter or signal '  param ' not found.']);
                        idxs(ip) = 0;
                    end
                    idxs(ip) = idx;
                end
            end

            params = params(logical(idxs));
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
            % BreachSet.CheckinDomain() Enforce parameters to adhere to their domains
            this.CheckinDomainParam();
            this.CheckinDomainTraj();
        end

        function CheckinDomainTraj(this)
            % BreachSet.CheckinDomainTraj()  Enforce signals to adhere to their domains

            if this.hasTraj()
                for itraj = 1:numel(this.P.traj)
                    for  i=1:this.P.DimX
                        if ~isequal(this.Domains(i).type, 'double')||~isempty(this.Domains(i).domain)
                            this.P.traj{itraj}.X(i,:) = this.Domains(i).checkin(this.P.traj{itraj}.X(i,:));
                        end
                    end
                end
            end
        end

        function CheckinDomainParam(this)
            pts = this.P.pts;
            if numel(this.Domains)< size(pts,1)
                this.Domains(size(pts,1)) = BreachDomain();
            end
            for i=this.P.DimX+1:size(pts,1)
                if ~isequal(this.Domains(i).type, 'double')||~isempty(this.Domains(i).domain)
                    pts(i,:) = this.Domains(i).checkin(pts(i,:));
                end
            end
            this.P.pts = pts;
        end

        %% Params
        function SetParam(this, params, values, is_spec_param)
            % BreachSet.SetParam(params, values [,  is_spec_param]) sets values to
            % parameters listed in params. If the set contains only one sample,
            % creates as many sample as there are values. If the set has
            % several samples and there is only one value, set this value to
            % all samples. Otherwise, returns an error.

            if nargin==2
                values = params;
                params = this.GetParamList();
            end


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
                elseif numel(values)==1
                    values = repmat(values, num_params,1);
                else
                    error('SetParam:wrong_arguments_size', 'Dimension mismatch between values and parameters.');
                end
            end

            num_pts  =  size(this.P.pts,2);
            num_values = size(values, 2);
            saved_traj = false;

            if ischar(is_spec_param)&&strcmp(is_spec_param, 'combine')&&...
                    ~(num_pts==1)&&~(num_values==1) % no need to combine (and mess order) in this case

                if this.hasTraj()
                    traj= this.P.traj;
                    saved_traj = true;
                    saved_Xf = this.P.Xf;
                end
                idx = N2Nn(2, [num_pts num_values]);
                old_pts = this.P.pts;
                this.P = Sselect(SPurge(this.P),1);
                this.P.pts = old_pts(:, idx(1,:));
                this.P.epsi= repmat(this.P.epsi,1, size(idx, 2));
                this.P.selected = zeros(1, size(idx, 2));
                if saved_traj
                    this.P.traj = traj;
                    this.P.Xf = saved_Xf;
                end
                this.P = SetParam(this.P, params, values(:, idx(2,:)));

            elseif ischar(is_spec_param)&&strcmp(is_spec_param, 'append')

                new_values = [GetParam(this.P, params) values];
                this.P.pts = [this.P.pts repmat(this.P.pts(:,end),1, size(values, 2))];
                this.P.epsi= [this.P.epsi repmat(this.P.epsi(:,end),1, size(values, 2))];
                this.P.selected = zeros(1, size(new_values, 2));
                this.P = SetParam(this.P, params, new_values);

            else  % legacy, i.e., not combine version
                if num_values==1 || num_values == num_pts
                    this.P = SetParam(this.P, params, values);
                elseif num_pts==1    % note in this case, we have to remove traces ( or see if maybe not, )
                    this.P.pts = repmat(this.P.pts,1, size(values, 2));
                    this.P.epsi= repmat(this.P.epsi,1, size(values, 2));
                    this.P.selected = zeros(1, size(values, 2));
                    this.P = SetParam(this.P, params, values);
                else
                    error('SetParam:wrong_arguments_size', 'Dimension mismatch between values and parameters.');
                end
            end

            this.ApplyParamGens(params);

        end

        function AddParam(this, params,values)
            % short for SetParam(..., 'append')

            if nargin == 2
                values = params;
                params = this.GetParamList;
            end
            this.SetParam(params, values, 'append');

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

        function SetDomainCfg(this, cfg)
            for ip = 1:numel(cfg.params)
                val = cfg.values{ip};
                if ischar(val)
                    val = str2num(val);
                end
                typ = cfg.types{ip};
                dom = cfg.domains{ip};
                if isempty(dom)
                    dom = [];
                elseif ischar(dom)
                    dom = str2num(dom); %#ok<ST2NM>
                elseif iscell(dom)
                    dom = cell2mat(dom);
                end
                this.SetParam(cfg.params{ip}, val);
                this.SetDomain(cfg.params{ip},typ,dom);
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
            if ~isempty(pg.domain)
                this.SetDomain(pg.params, pg.domain);
            else
                domain = repmat(BreachDomain, 1, numel(pg.params));
                this.SetDomain(pg.params, domain);
            end

            
            % this.SetParam(pg.params_out, nan, true); % why ?
            params_nan = setdiff(pg.params_out, pg.params);
            if ~isempty(params_nan)
                this.SetParam(params_nan, nan, true);
            end
            
            % update domain of output parameters
            if ~isempty(pg.domain_out)
                for ip =1:numel(pg.params_out)
                    this.SetDomain(pg.params_out{ip}, pg.domain_out(ip));
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
                if ~isempty(ig)
                    this.P = Preset_traj_ref(this.P);
                end
            end
        end

        function values = GetParam(this, params, ip)
            if nargin==1||isempty(params)
                params= this.GetParamList;
            end

            values = GetParam(this.P,params);

            if ~isempty(values)&&nargin>=3
                values = values(:, ip);
            end
        end

        function [values, i1, i2] = GetUParam(this,params,argu)
            % [values, i1, i2] = GetUParam(this, params,argu) returns unique values of parameters
            % i1 are indices of unique values, i2 is a cell of (unique) indices
            % for where the unique values were found. argu is an argument for
            % unique function call (for 'stable' essentially)
            %

            if nargin<3
                argu = 0;
            end

            all_values = this.GetParam(params);
            if argu
                [values, i1, ib] = unique(all_values','rows', argu);
            else
                [values, i1, ib] = unique(all_values','rows');
            end

            values = values';
            if nargout >= 1
                i1 = i1';
            end

            if nargout >=2
                ib = ib';
                i2 = cell(1,numel(i1));
                for ind = 1:numel(i1)
                    i2{ind} = find(ib==ind);
                end
            end
        end

        function RemoveDuplicateParams(this)
            params = this.GetVariables();
            [~, i1, ~] = GetUParam(this,params,'stable');
            this.P = Sselect(this.P, i1);
        end

        function RemoveParams(this,params, values)
            val_before = this.GetParam(params);
            [~, idx_after] = setdiff(round(val_before,14)',round(values,14)',"rows");
            if ~isempty(idx_after)
                this.P = Sselect(this.P, idx_after);
            else
                this.ResetParamSet();
                this.P.pts = this.P.pts(:, []);
                this.P.epsi  = this.P.epsi(:,[]);                
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
                    this.Domains(i_params(ip)).domain = ranges(ip,:);
                    %warning('SetParamRanges:enum_or_bool', 'Use SetDomain
                    %for enum or bool types.' ); % Maybe should keep the
                    %warning
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

        % Get the number of param samples - -1 means P is empty
        function nb_pts = GetNbParamVectors(this)
            if isempty(this.P)
                nb_pts = -1;
            else
                nb_pts= size(this.P.pts,2);
            end
        end
        % Same
        function nb_pts= GetNumSamples(this)
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
                [params, i_diff] = setdiff(params, req_params, 'stable');
                ipr = ipr(i_diff);
            end
        end

        function [params, ipr] = GetReqVariables(this)
            [params, ipr] = GetVariables(this);
            if ~isempty(params)
                req_params = this.GetPropParamList();
                [params, i_intersect] = intersect(params, req_params, 'stable');
                ipr = ipr(i_intersect);
            end
        end

        function [ params, ipr]  = GetBoundedDomains(this)
            % GetNonEmptyDomains
            ipr = cellfun(@(c)(~isempty(c)), {this.Domains.domain});
            ipr = find(ipr);
            params =   this.P.ParamList(ipr);
        end

        function SetEmptyDomain(this, params)

            if ~iscell(params)
                params = {params};
            end
            for idx_param = 1:numel(params)
                idx_param_in_this = FindParam(this.P, params{idx_param});
                this.Domains(idx_param_in_this).domain = [];
            end

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
                            sig1 = varargin{1}{is};
                            sig2 = varargin{2}{is};
                            add_sigs(sig1, sig2);
                        end
                    else
                        sig1 = varargin{1};
                        sig2 = varargin{2};
                        add_sigs(sig1, sig2);
                    end
                otherwise
                    for is = 1:numel(varargin)/2
                        try
                            sig1 = varargin{2*is-1};
                            sig2 = varargin{2*is};
                            add_sigs(sig1, sig2);
                        catch
                            error('SetSignalMap:wrong_arg', arg_err_msg);
                        end
                    end
            end

            if this.verbose >= 2
                this.PrintSigMap();
            end

            function add_sigs(sig1, sig2)
                if ischar(sig1)&&ischar(sig2)
                    if ~strcmp(sig1,sig2)
                        this.sigMap(sig1) = sig2;
                        this.sigMapInv(sig2) = sig1;

                        % get current aliases for sig1 and sig2
                        if this.AliasMap.isKey(sig1)
                            SIG1 = [{sig1} this.AliasMap(sig1)];
                        else
                            SIG1 = {sig1};
                        end

                        if this.AliasMap.isKey(sig2)
                            SIG2 = [{sig2} this.AliasMap(sig2)];
                        else
                            SIG2 = {sig2};
                        end

                        % add aliases of each other
                        if ~this.AliasMap.isKey(sig1)
                            this.AliasMap(sig1) = SIG2;
                        else
                            this.AliasMap(sig1) = unique([this.AliasMap(sig1) SIG2], 'stable');
                        end
                        if ~this.AliasMap.isKey(sig2)
                            this.AliasMap(sig2) = SIG1;
                        else
                            this.AliasMap(sig2) = unique([this.AliasMap(sig2) SIG1], 'stable');
                        end
                    end
                else
                    error('SetSignalMap:wrong_arg', arg_err_msg);
                end
            end

        end

        function ResetSigMap(this)
            this.sigMap = containers.Map();
            this.sigMapInv = containers.Map();
            this.AliasMap =containers.Map();
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
                traj = this.P.traj{itrajs(i_traj)};
                Xi = traj.X;
                if (~exist('t','var'))
                    X{i_traj} = Xi(signals_idx,:);
                else
                    X{i_traj} = interp1(traj.time, Xi(signals_idx,:)',t)';
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

            idx = zeros(1, numel(signals));
            ifound = idx;
            for isig = 1:numel(signals)
                sig = signals{isig};
                [idx(isig), ifound(isig)] = FindParam(this.P, sig);
                if this.AliasMap.isKey(sig)
                    aliases_sig = this.AliasMap(sig);
                    for ais = 1:numel(aliases_sig)
                        [idx_s, ifound_s] = FindParam(this.P, aliases_sig{ais});
                        if ifound_s
                            ifound(isig)=true;
                            idx(isig)=idx_s;
                        end
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

        function h= PlotSurf(this)
            
            
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
            elseif isnumeric(params)
                switch nargin
                    case 2 % Sample randomly
                    case 3
                        method = num_samples;
                    case 4
                        opt_multi = method;
                        method = num_samples;
                    case 5
                        max_num_samples = opt_multi;
                        opt_multi = method;
                        method = num_samples;
                end
                num_samples =params;
                params = this.GetVariables();
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

            if isempty(domains)
                error('sample:empty_domain', 'Domain to sample is empty or undefined.');
            end


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
                this.P = Pimport_traj(this.P, P0);
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
            this.ApplyParamGens();
        end

        function PseudoRandomSample(this, nb_sample)
            % Pseudo-random sampling

            this.ResetParamSet();

            newP = TestronRefine(this.P, nb_sample);

            this.P = newP;
            this.CheckinDomainParam();
            this.ApplyParamGens();
        end

        function Bm = MorrisSample(this, vars, ranges, num_path, size_grid, seed)

            if isempty(vars)
                vars = this.GetVariables();
            end

            if isempty(ranges)
                ranges = this.GetParamRanges(vars);
            end

            if isequal(size(ranges), [2 1])
                ranges = [1 2];
            end

            if isequal(size(ranges), [1 2])
                ranges = repmat(ranges, numel(vars),1);
            end


            Bm = BreachSet(vars);
            Bm.SetParamRanges(vars, ranges);
            Bm.P.epsi = ((ranges(:,2)+ranges(:,1))/2)'; % legacy stuff
            Pr = pRefine(Bm.P, size_grid,num_path,seed);
            X0 = Pr.pts;
            this.ResetParamSet;
            this.SetParam(vars, X0)
            this.P.opt_morris =struct('num_path',num_path,'size_grid',size_grid,'rand_seed',1);
            this.P.D = Pr.D;
            if nargout>=1
                Bm.P = Pr;
            end
        end

        %% Concatenation, ExtractSubset - needs some additional compatibility checks...
        function Concat(this, other, fast)
            if nargin<=2
                fast = false;
            end
            this.P = SConcat(this.P, other.P, fast);
        end

        function other  = ExtractSubset(this, idx)
            other = this.copy();
            other.P = Sselect(this.P, idx);
        end

        function SavedTrajectorySample(this, paramValues)
            % TESTRON
            % Fixes the "sampling" for a stored trajectory

            this.ResetParamSet();

            newP = SavedTrajectoryRefine(this.P, paramValues);

            this.P = newP;
            this.CheckinDomainParam();
        end



        %% Plot parameters
        function ax= PlotPts(this, params, ax, varargin)
            % Simple plot function.
            if ~exist('params','var')
                params = this.GetParamList();
            end

            if ~exist('ax','var')
                ax = gca;
            end

            if isempty(varargin)
                varargin = {'bx'};
            end

            switch numel(params)
                case 1
                    values = this.GetParam(params);
                    axes(ax);
                    plot(values, 0*values, varargin{:});

                case 2
                    values = this.GetParam(params);
                    axes(ax);
                    plot(values(1,:), values(2,:), varargin{:});
                otherwise
                    params =params(1:3);
                    values = this.GetParam(params);
                    axes(ax);
                    plot3(values(1,:), values(2,:), values(3,:), varargin{:});
            end
        end

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

            function [txt] = myupdatefcn(~,event_obj)

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

        function ax = PlotBoxPts(this, varargin)
            % legacy plot box function. Rescucitated for grids returned by
            % GridFilter and GridFilterSignals
            ax = gca;
            SplotBoxPts(this.P, varargin{:});
            
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

        function cfg=  GetDomainCfg(this)
            params = this.GetParamList();
            cfg.params = params;
            for ip = 1:numel(params)
                dom = this.GetDomain(params{ip});
                cfg.types{ip} = dom.type;
                cfg.values{ip} = this.GetParam(params{ip},1);
                if strcmp(dom.type, 'enum')
                    cfg.domains{ip}=dom.enum;
                else
                    cfg.domains{ip}=dom.domain;
                end
            end
        end

        %% Coverage

        function [cover_opts] = SetCoverageOptions(this,cover_opts, varargin)
            % Create a structure objects for grid filtering and coverage computation options

            other_opts.ParamsMode = 'auto';          % 'auto' detects variable parameters and sets ranges manually if specified or
            other_opts.SignalsMode = 'auto';          % 'off' means don't compute
            other_opts.SignalsFreqsMode = 'off';    % freq domain off by default for now

            other_opts.IncludeParams = 'all';
            other_opts.ExcludeParams = {};

            other_opts.IncludeSignals = 'all';
            other_opts.ExcludeSignals = {};

            other_opts.GridsMode = 'auto';  % 'auto' means try to guess if it is a resolution or a size
            % other choices are
            % 'size_cells', 'num_cells'


            if nargin>1&&ischar(cover_opts)
                other_opts = varargin2struct_breach(other_opts,cover_opts, varargin{:});
                cover_opts = [];
            else
                other_opts = varargin2struct_breach(other_opts, varargin{:});
            end

            if ~exist('cover_opts', 'var')||isempty(cover_opts)
                cover_opts =struct();
                grid_spec = 10;
            elseif isnumeric(cover_opts)
                grid_spec = cover_opts;
                cover_opts = struct;
            else
                grid_spec = 10; % default grid spec for when it is not specified by cover_opts already
            end

            %% parameters coverage
            if ~isequal(other_opts.ParamsMode, 'off' )

                % if parameters are specified in opts, take those,
                % otherwise take all variables
                if isfield(cover_opts, 'params')
                    params = fieldnames(cover_opts.params)';
                elseif isequal(other_opts.IncludeParams,'all')
                    params = this.GetVariables();
                else
                    params = this.expand_param_name(other_opts.IncludeParams);
                end

                if ~isempty(other_opts.ExcludeParams)
                    exclude_list = this.expand_param_name(other_opts.ExcludeParams);
                    params = setdiff(params,exclude_list, 'stable');
                end
                
                [idx_params, found] = FindParam(this.P,params);
                if ~isempty(params)
                    if ~isfield(cover_opts,'params')
                        cover_opts.params = struct();
                    end

                    for ip = 1:numel(params)
                        this_param = params{ip};
                        if ~found(ip)
                            warning('Parameter %s not found, ignored', this_param);
                            continue
                        end
                        idx_pts = idx_params(ip);
                        this_domain = this.Domains(idx_pts);   % existing breach_domain for this param

                        % add the parameter if not in opts already
                        if ~isfield(cover_opts.params, this_param)
                            cover_opts.params.(this_param)  =struct();
                        end

                        % checks its range
                        if ~isfield(cover_opts.params.(this_param),'range')
                            if ~isempty(this_domain.domain)
                                range = this_domain.domain;
                            else
                                Rmin = min(this.P.pts(idx_pts,:));  % could be more efficient with bool and enum but whatever
                                Rmax = max(this.P.pts(idx_pts,:));
                                range = [Rmin, Rmax];
                            end
                        else
                            range = cover_opts.params.(this_param).range;
                        end

                        % checks grid specifications - we get the grid, then
                        % makes sure it intersects with the range (if it was
                        % specified)
                        if isfield(cover_opts.params.(this_param), 'grid')
                            grid =cover_opts.params.(this_param).grid;
                        elseif iscell(grid_spec)
                            try
                                grid = grid_spec{ip};
                            catch
                                warning('Problem with grid_spec dimensions, using grid_spec(1)');
                                grid = grid_spec{1};
                            end
                        elseif isscalar(grid_spec)
                            grid = grid_spec;
                        else
                            try
                                grid = grid_spec(ip);
                            catch
                                warning('Problem with grid_spec dimensions, trying grid_spec(1)');
                                grid = grid_spec(1);
                            end
                        end
                        cover_opts.params.(this_param).grid = grid;
                        % if enum is not included in range, intersect it
                        % this is where we checks consistency with original
                        % type if defined
                        switch this_domain.type
                            case 'bool'
                                enum = [0 1];
                                range = [0 1];
                                this_domain.domain = [0 1]; % not automatic I believe
                                eps = enum*0+.1;
                            case 'int' % here we consider the grid should be all integers between min and max, regardless of grid spec
                                this_domain.domain = range;
                                this_domain.enum = range(1):range(2);
                                enum = this_domain.enum;
                                eps = enum*0+.1;
                            case 'enum' % we exclude the enum outside the range
                                this_domain.domain = range;
                                if min(this_domain.enum) < range(1)
                                    this_domain.enum = this_domain.enum(this_domain.enum>= range(1));
                                end
                                if max(this_domain.enum) > range(2)
                                    this_domain.enum = this_domain.enum(this_domain.enum <= range(2));
                                end
                                enum = this_domain.enum;
                                eps = min(diff(enum))/10;
                                eps = repmat(eps, 1, numel(enum));
                            case 'double'
                                % translates grid spec (ex delta) into enum
                                switch other_opts.GridsMode
                                    case 'auto'
                                        if isscalar(grid)&&grid==floor(grid)&&grid>1
                                            eps = (range(2)-range(1))/(grid*2);
                                            enum = linspace(range(1)+eps,range(2)-eps, grid); % delta samples in dim ip
                                            eps = repmat(eps, 1, numel(enum));
                                        elseif isscalar(grid)&&grid>0
                                            eps = min([grid/2, (range(2)-range(1))/2]);
                                            num_cell = floor((range(2)-range(1))/grid); 
                                            start_enum = range(1)+((range(2)-range(1)) -num_cell*grid)/2 ;  % ensure the grid is centered wrt range, inside range and cells cover range
                                            enum = start_enum:grid:range(2);               % grid is resolution
                                            eps = repmat(eps, 1, numel(enum));
                                        else
                                            % we need to adjust (center)
                                            eps = diff(grid)/2;
                                            enum = grid(2:end)-eps;
                                        end
                                    case 'num_cells'
                                        eps = (range(2)-range(1))/(grid*2);
                                        enum = linspace(range(1)+eps,range(2)-eps, grid); % delta samples in dim ip
                                        eps = repmat(eps, 1, numel(enum));

                                    case 'size_cells'
                                        eps = grid/2;
                                        num_cell = floor((range(2)-range(1))/grid); 
                                        start_enum = range(1)+((range(2)-range(1)) -num_cell*grid)/2 ;  % ensure the grid is centered wrt range, inside range and cells cover range
                                        enum = start_enum:grid:range(2);               % grid is resolution
                                        eps = repmat(eps, 1, numel(enum));

                                    otherwise
                                        error('Wrong mode ''%s'' for GridsMode, should be ''auto'',''num_cells'' or ''size_cells''', other_opts.GridsMode);
                                end

                                this_domain = BreachDomain('enum', range, enum);
                        end

                        % feed back enum as grid -- not. ( problematic with
                        % custom grid) 
                        cover_opts.params.(this_param).range = range;
                        %cover_opts.params.(this_param).grid = enum;
                        cover_opts.params.(this_param).eps = eps;
                        cover_opts.params.(this_param).breach_domain = this_domain;

                    end
                    for size_proj = 1:min(3, numel(params)) % default max projection size is 3
                        cover_opts.projections{size_proj} = get_projections(params,size_proj);
                    end
                end
            else
                if isfield(cover_opts,'params')
                    cover_opts = rmfield(cover_opts,'params');
                end
            end


            %%  Signals
            if ~isequal(other_opts.SignalsMode, 'off' )&&this.hasTraj()


                % if signals are specified in opts, take those,
                % otherwise take all
                
                if isfield(cover_opts, 'signals')
                    signals = fieldnames(cover_opts.signals)';
                elseif isequal(other_opts.IncludeSignals,'all')
                    signals = this.GetSignalList();
                else
                    signals = this.expand_signal_name(other_opts.IncludeSignals);
                end

                if ~isempty(other_opts.ExcludeSignals)
                    exclude_list = this.expand_signal_name(other_opts.ExcludeSignals);
                    signals = setdiff(signals,exclude_list, 'stable');
                end

                % time range
                if ~isfield(cover_opts,'signals')
                    cover_opts.signals = struct();
                end
                if ~isfield(cover_opts.signals, 'time')
                    cover_opts.signals.time = struct();
                end

                % default time range
                time = this.GetTime();
                time_range= [time(1) time(end)];
                if ~isfield(cover_opts.signals.time,'range')
                    cover_opts.signals.time.range = time_range;
                else
                    range = cover_opts.signals.time.range;
                    if size(range,2)~=2
                        warning('Invalid time range, using default.');
                        cover_opts.signals.time.range = time_range;
                    end
                end

                %% time grid
                if isfield(cover_opts.signals.time, 'grid')
                    time_grid = cover_opts.signals.time.grid;
                else
                    time_grid = 10;
                end

                if isscalar(time_grid)&&time_grid==floor(time_grid)&&time_grid>1
                    eps = (time_range(2)-time_range(1))/(time_grid*2);
                    enum = linspace(time_range(1)+eps,time_range(2)-eps, time_grid); % delta samples in dim ip
                    eps = repmat(eps, 1, numel(enum));
                elseif isscalar(time_grid)&&time_grid>0
                    eps = time_grid/2;
                    num_cell = floor((range(2)-range(1))/time_grid); 
                    start_enum = range(1)+((range(2)-range(1)) -num_cell*time_grid)/2 ;  % ensure the grid is centered wrt range, inside range and cells cover range
                    enum = start_enum:time_grid:range(2);                                                 % grid is resolution
                    eps = repmat(eps, 1, numel(enum));
                else
                    eps = diff(time_grid)/2;
                    enum = time_grid;
                end

                cover_opts.signals.time.range = time_range;
                cover_opts.signals.time.grid = enum;
                cover_opts.signals.time.eps = eps;
                cover_opts.signals.time.breach_domain =  BreachDomain('enum', time_range, enum);

                %% signals range and grid
                [~, found_signals] = FindParam(this.P, signals);
                for ip = 1:numel(signals)
                    this_signal = signals{ip};
                    if ~found_signals(ip)
                        warning('Signal %s not found, ignore and continue (at your own risk).', this_signal);
                    end

                    this_domain = this.GetDomain(this_signal);   % existing breach_domain for this signal

                    % add the signal if not in opts already
                    if ~isfield(cover_opts.signals, this_signal)
                        cover_opts.signals.(this_signal)  =struct();
                    end

                    % checks its range
                    if ~isfield(cover_opts.signals.(this_signal),'range')
                        if ~isempty(this_domain.domain)
                            range = this_domain.domain;
                        else
                            sig_values = this.GetSignalValues(this_signal);
                            if ~iscell(sig_values)
                                sig_values = {sig_values};
                            end
                            Rmin = inf;
                            Rmax  =-inf;
                            for itraj = 1:numel(sig_values)
                                Rmin = min(Rmin, min(sig_values{itraj}));
                                Rmax = max(Rmax, max(sig_values{itraj}));
                            end
                            range = [Rmin, Rmax];
                        end
                    else
                        range = cover_opts.signals.(this_signal).range;
                    end

                    % checks grid specifications - we get the grid, then
                    % makes sure it intersects with the range (if it was
                    % specified)
                    if isfield(cover_opts.signals.(this_signal), 'grid')
                        grid =cover_opts.signals.(this_signal).grid;
                    elseif iscell(grid_spec)
                        try
                            grid = grid_spec{ip};
                        catch
                            warning('Problem with grid_spec dimensions for %s, trying grid_spec{1}',this_signal);
                            grid = grid_spec{1};
                        end
                    elseif isscalar(grid_spec)
                        grid = grid_spec;
                    else
                        try
                            grid = grid_spec(ip);
                        catch
                            warning('Problem with grid_spec dimensions for %s, trying grid_spec{1}',this_signal);
                            grid = grid_spec(1);
                        end
                    end


                    % if enum is not included in range, intersect it
                    % this is where we checks consistency with original
                    % type if defined
                    switch this_domain.type
                        case 'bool'
                            enum = [0 1];
                            range = [0 1];
                            this_domain.domain = [0 1]; % not automatic I believe
                        case 'int' % here we consider the grid should be all integers between min and max, regardless of grid spec
                            this_domain.domain = range;
                            this_domain.enum = range(1):range(2);
                            enum = this_domain.enum;
                            eps = enum*0+.1;
                        case 'enum' % we exclude the enum outside the range
                            this_domain.domain = range;
                            if min(this_domain.enum) < range(1)
                                this_domain.enum = this_domain.enum(this_domain.enum>= range(1));
                            end
                            if max(this_domain.enum) > range(2)
                                this_domain.enum = this_domain.enum(this_domain.enum <= range(2));
                            end
                            enum = this_domain.enum;
                            eps = min(diff(enum))/10;
                            eps = repmat(eps, 1, numel(enum));
                        case 'double'
                            % translates grid spec (ex delta) into enum
                            switch other_opts.GridsMode
                                case 'auto'
                                    if isscalar(grid)&&grid==floor(grid)&&grid>1
                                        eps = (range(2)-range(1))/(grid*2);
                                        enum = linspace(range(1)+eps,range(2)-eps, grid); % delta samples in dim ip
                                        eps = repmat(eps, 1, numel(enum));
                                    elseif isscalar(grid)&&grid>0
                                        eps = grid/2;
                                        num_cell = floor((range(2)-range(1))/grid); 
                                        start_enum = range(1)+((range(2)-range(1)) -num_cell*grid)/2 ;  % ensure the grid is centered wrt range, inside range and cells cover range                                           
                                        enum = start_enum:grid:range(2);               % grid is resolution
                                        eps = repmat(eps, 1, numel(enum));
                                    else
                                        % we need to adjust (center)
                                        eps = diff(grid)/2;
                                        enum = grid(2:end)-eps;
                                    end
                                case 'num_cells'
                                    eps = (range(2)-range(1))/(grid*2);
                                    enum = linspace(range(1)+eps,range(2)-eps, grid); % delta samples in dim ip
                                    eps = repmat(eps, 1, numel(enum));

                                case 'size_cells'
                                    eps = grid/2;
                                    num_cell = floor((range(2)-range(1))/grid); 
                                    start_enum = range(1)+((range(2)-range(1)) -num_cell*grid)/2 ;  % ensure the grid is centered wrt range, inside range and cells cover range                                          
                                    enum = start_enum:grid:range(2);               % grid is resolution                   
                                    eps = repmat(eps, 1, numel(enum));

                                otherwise
                                    error('Wrong mode ''%s'' for GridsMode, should be ''auto'',''num_cells'' or ''size_cells''', other_opts.GridsMode);
                            end
                            this_domain = BreachDomain('enum', range, enum);
                    end

                    % feed back enum as grid (or not ? )
                    cover_opts.signals.(this_signal).range = range;
                    cover_opts.signals.(this_signal).grid = enum;
                    cover_opts.signals.(this_signal).eps =eps;
                    cover_opts.signals.(this_signal).breach_domain = this_domain;

                end

                % Projections
                for size_proj = 1:min(3, numel(signals)) % default max projection size is 3
                    cover_opts.signals_projections{size_proj} = get_projections(signals,size_proj);
                end

                for size_proj = 1:min(3, numel(signals)) % time projections, append time before the above
                    cover_opts.time_signals_projections{size_proj} = cellfun( @(c) (['time' c ]), cover_opts.signals_projections{size_proj}  , 'UniformOutput', false);
                end
            else
                if isfield(cover_opts,'signals')
                    cover_opts = rmfield(cover_opts,'signals');
                end
            end

            if ~isequal(other_opts.SignalsFreqsMode, 'off' )&&this.hasTraj()
                % freqs ranges and grid
                if ~isfield(cover_opts, 'signals_freqs')
                    cover_opts.signals_freqs = struct();
                end

                % get target frequency range from time
                traj = this.P.traj{1}; % assumes uniform traces....
                Ts = min(diff(traj.time));
                time = 0:Ts:traj.time(end);
                if rem(numel(time),2)
                    time = time(1:end-1); % ensures pair number of samples for FFT
                end
                Fs = 1/Ts;
                %L = numel(time);
                fs_max = min(100, Fs/2); % Max a hundred hertz unless time step is too large for that
                freq_range = [0 fs_max];

                signals=this.GetSignalList();
                for ip = 1:numel(signals)
                    this_signal = [signals{ip} '_freqs'];

                    if ~isfield(cover_opts.signals_freqs,this_signal)
                        cover_opts.signals_freqs.(this_signal) = struct();
                    end

                    % range: default unless specified otherwise
                    if ~isfield(cover_opts.signals_freqs.(this_signal),'range')
                        cover_opts.signals_freqs.(this_signal).range = freq_range;
                    else
                        freq_range = cover_opts.signals_freqs.(this_signal).range;
                    end

                    % grid spec - as for double domains, and defaults to 10
                    if isfield(cover_opts.signals_freqs.(this_signal), 'grid')
                        freq_grid = cover_opts.signals_freqs.(this_signal).grid;
                    elseif iscell(grid_spec)
                        try
                            freq_grid = grid_spec{ip};
                        catch
                            warning('Problem with grid_spec dimensions for %s, trying grid_spec{1}',this_signal);
                            freq_grid = grid_spec{1};
                        end
                    elseif isscalar(grid_spec)
                        freq_grid = grid_spec;
                    else
                        try
                            freq_grid = grid_spec(ip);
                        catch
                            warning('Problem with grid_spec dimensions for %s, trying grid_spec{1}',this_signal);
                            freq_grid = grid_spec(1);
                        end
                    end


                    switch other_opts.GridsMode
                        case 'auto'

                            if isscalar(freq_grid)&&freq_grid==floor(freq_grid)&&freq_grid>1
                                eps = (freq_range(2)-freq_range(1))/(freq_grid*2);
                                enum = linspace(freq_range(1)+eps,freq_range(2)-eps, freq_grid); % delta samples in dim ip
                                eps = repmat(eps, 1, numel(enum));
                            elseif isscalar(freq_grid)&&freq_grid>0                                
                                eps = freq_grid/2;                                
                                num_cell = floor((range(2)-range(1))/freq_grid);
                                start_enum = range(1)+((range(2)-range(1)) -num_cell*freq_grid)/2 ;  % ensure the grid is centered wrt range, inside range and cells cover range
                                enum = start_enum:freq_grid:range(2);               % grid is resolution
                                eps = repmat(eps, 1, numel(enum));
                            else
                                eps = diff(freq_grid)/2;
                                enum = freq_grid;
                            end
                        case 'num_cells'
                            eps = (freq_range(2)-freq_range(1))/(freq_grid*2);
                            enum = linspace(freq_range(1)+eps,freq_range(2)-eps, freq_grid); % delta samples in dim ip
                            eps = repmat(eps, 1, numel(enum));

                        case 'size_cells'
                            eps = freq_grid/2;
                            num_cell = floor((range(2)-range(1))/freq_grid); 
                            start_enum = range(1)+((range(2)-range(1)) -num_cell*freq_grid)/2 ;  % ensure the grid is centered wrt range, inside range and cells cover range                                           
                            enum = start_enum:freq_grid:range(2);               % grid is resolution
                            eps = repmat(eps, 1, numel(enum));

                        otherwise
                            error('Wrong mode ''%s'' for GridsMode, should be ''auto'',''num_cells'' or ''size_cells''', other_opts.GridsMode);
                    end

                    cover_opts.signals_freqs.(this_signal).range = freq_range;
                    cover_opts.signals_freqs.(this_signal).grid = enum;
                    cover_opts.signals_freqs.(this_signal).eps = eps;
                    cover_opts.signals_freqs.(this_signal).breach_domain = BreachDomain('enum', freq_range,enum);

                end

                signals_for_freqs_coverage = fieldnames(cover_opts.signals_freqs)';

                for size_proj = 1:min(3, numel(signals_for_freqs_coverage))
                    cover_opts.signals_freqs_projections{size_proj} = get_projections(signals_for_freqs_coverage,size_proj);
                end

            end

            this.CoverOpts = cover_opts;

            %% Auxiliary functions
            function projs = get_projections(params, n)

                if n==1
                    projs= cell(1, numel(params));
                    for idx_p= 1:numel(params)
                        projs{idx_p} =params(idx_p);
                    end
                else
                    set_params_idx = nchoosek(1:numel(params),n);
                    size_set = size(set_params_idx,1);
                    projs = cell(1,size_set);
                    for idx_set= 1:size_set  % for each combination
                        this_set_params_idx = set_params_idx(idx_set,:);
                        this_proj = cell(1, n);
                        for idx_p = 1:n
                            this_proj{idx_p} = params{this_set_params_idx(idx_p)};
                        end
                        projs{idx_set} = this_proj; % parameters in this combination
                    end
                end
            end


        end

        function [cov, Bg, Bs, Bgf] = GetCoverage(this, opts)
            % Returns coverage structure, and occupancy grid for
            % parameters and signals

            % Checks options. If empty or not given, checks this.CoverOpts
            opts_is_same_as_this = false;
            if nargin<=1||isempty(opts)
                if ~isempty(this.CoverOpts)
                    opts = this.CoverOpts;
                    opts_is_same_as_this = true;
                else
                    opts = this.SetCoverageOptions();
                end
            end

            if ~opts_is_same_as_this % not sure yet
                if isequal(opts,this.CoverOpts)
                    opts_is_same_as_this = true; % in fact it is
                end
            end

            this.CoverOpts = opts; % ok, will go with that.

            % checks if we have previous results, and they are consistent
            if  opts_is_same_as_this&&~isempty(this.CoverRes)

                pts_prev = this.P.pts(:, 1:this.CoverRes.num_samples);
                hash_prev = DataHash(pts_prev);
                if ~isequal(hash_prev, this.CoverRes.hash_samples)
                    this.CoverRes=[];  % not same as before, start over
                end
            end

            % previous results ?
            if ~isempty(this.CoverRes)
                cov = this.CoverRes;
            else
                cov  = struct();
                Bg = [];
                Bs = [];
                Bgf = [];
            end

            if isfield(opts,'params')
                if ~isequal(cov, struct())% we got previous results
                    try
                        Bg = this.CoverRes.params.Bg;
                    catch
                        warning('Previous grid set absent ?'); % should not happen, for debug
                        Bg= [];
                    end
                end
                
                Bg = this.GridFilter([], opts, Bg);
                params= setdiff(Bg.GetVariables(),{'traj__ref__','count','idx'},'stable');
                
                for iproj_set  =1:numel(opts.projections) % loop over dimensions of projections
                    cov_dim = 0;
                    proj_set = opts.projections{iproj_set};

                    num_proj = numel(proj_set);
                    for iproj = 1:num_proj
                        this_proj= proj_set{iproj};
                        this_cov = compute_coverage(this_proj);
                        cov_dim  = cov_dim + this_cov;
                    end
                    cov_dim = cov_dim/num_proj;
                    cov.params.(['cov' num2str(iproj_set) 'd']) = cov_dim;
                    cov.params.Bg = Bg;
                end

            end

            % signals
            if isfield(opts,'signals')

                %% Todo incremental coverage from previous result
                Bs_all  = this.sig2params();
                opts_sigs.params = opts.signals;
                opts_sigs.projections = opts.signals_projections;
                num_sig_proj = numel( opts.signals_projections);

                timed_projs = opts.time_signals_projections;
                num_time_sig_proj = numel( opts.time_signals_projections);

                opts_sigs.projections = [opts_sigs.projections timed_projs];
                [cov_signals, Bs]  = Bs_all.GetCoverage(opts_sigs);
                cov.signals = cov_signals.params;

                % rename aggregated measures
                for ip  = 1:num_sig_proj
                    old_name = ['cov' num2str(ip) 'd'];
                    name = ['cov' num2str(ip) 'd_sig'] ;
                    cov.signals.(name) = cov.signals.(old_name);
                    cov.signals = rmfield(cov.signals,old_name);
                end

                for ip  = 1:num_time_sig_proj
                    old_name = ['cov' num2str(ip+num_sig_proj) 'd'];
                    name = ['cov' num2str(ip) 'd_time_sig'] ;
                    cov.signals.(name) = cov.signals.(old_name);
                    cov.signals = rmfield(cov.signals,old_name);
                end
                cov.signals.Bg = Bs;
                cov.signals.Bs_all = Bs_all;
            end

            if isfield(opts,'signals_freqs')
                sigs_freqs = fieldnames(opts.signals_freqs)';
                Bf = BreachSet(sigs_freqs);
                for i=1:numel(sigs_freqs)
                    this_freq_sig = sigs_freqs{i} ;
                    Bf.SetParamRanges(this_freq_sig, opts.signals_freqs.(this_freq_sig).range);
                end
                sigs_names_for_freq_coverage =  cellfun(@(c)(c(1:end-6)), sigs_freqs, 'UniformOutput',false);
                sigs_values_for_freq_coverage = this.GetSignalValues(sigs_names_for_freq_coverage);
                if ~iscell(sigs_values_for_freq_coverage)
                    sigs_values_for_freq_coverage={sigs_values_for_freq_coverage};
                end

                for i=1:numel(sigs_freqs)
                    this_freq_sig = sigs_freqs{i};
                    cov_opts_for_freq_coverage.params.(this_freq_sig) = opts.signals_freqs.(this_freq_sig);
                    cov_opts_for_freq_coverage.projections = opts.signals_freqs_projections;
                end

                %  Apply FFT
                for itraj = 1:numel(sigs_values_for_freq_coverage)
                    Bf_itraj = BreachSet(sigs_freqs);
                    time = this.P.traj{itraj}.time;
                    Ts = time(2)-time(1);
                    if any(diff(time)-Ts>1e-10)
                        error('Yup, irregular sampling time, time to implement this time interpolation again.');
                    end
                    %time = 0:Ts:traj.time(end);
                    %if rem(numel(time),2)
                    %    time = time(1:end-1); % ensures pair number of samples for FFT
                    %end
                    Fs = 1/Ts;
                    L = numel(time);
                    X = sigs_values_for_freq_coverage{itraj}';
                    all_freqs = Fs*(0:(L/2))/L; % data resolution

                    dim_sigs =numel(sigs_freqs);
                    if size(X,2)~=dim_sigs
                        X=X';
                    end
                    for isig =1:dim_sigs
                        this_freq_sig = sigs_freqs{isig} ;
                        Y = fft(X(:,isig));
                        P2 = abs(Y/L);
                        X_hat = P2(1:floor(L/2)+1);
                        X_hat(2:end-1) =  2*X_hat(2:end-1);
                        X_hat_threshold = 0.1;

                        %MaxX = max(X_hat);
                        %X_hat_threshold = MaxX/10; % kind of arbitrary...
                        % TODO: options for threshold, or better default
                        % choice...

                        freqs = all_freqs(X_hat>X_hat_threshold);
                        DEBUG = 0;
                        if DEBUG
                            figure;
                            subplot(2,1,1)
                            plot(X(:,isig));
                            subplot(2,1,2)
                            plot(all_freqs, X_hat)
                            hold on
                            plot(all_freqs, 0*X_hat+X_hat_threshold);
                            plot(freqs, X_hat(X_hat>X_hat_threshold),'x-r');
                            pause;
                            close;

                        end
                        % collect frequencies for this trace and signal,
                        % filter them before combining
                        Bf_itraj_isig = BreachSet(this_freq_sig);
                        Bf_itraj_isig.SetParamRanges(this_freq_sig,opts.signals_freqs.(this_freq_sig).range);
                        this_opts.params = struct();
                        this_opts.params.(this_freq_sig)= opts.signals_freqs.(this_freq_sig);

                        Bf_itraj_isig.SetParam(this_freq_sig,  freqs);
                        Bg_itraj_isig =Bf_itraj_isig.GridFilter([], this_opts);
                        filtered_freqs = Bg_itraj_isig.GetParam(this_freq_sig);

                        if isig==1
                            Bf_itraj.SetParam(this_freq_sig, filtered_freqs);
                        else
                            Bf_itraj.SetParam(this_freq_sig,filtered_freqs, 'combine')
                        end

                    end

                    if itraj==1 % first trace
                        Bf=Bf_itraj;
                    else
                        Bf.Concat(Bf_itraj,1) % concat freqs for subsequent traces
                    end
                end

                % apply coverage on BF and copy to freq_signals stuff
                [cov_for_freq_coverage, Bgf] = Bf.GetCoverage(cov_opts_for_freq_coverage);
                cov.signals_freqs = cov_for_freq_coverage.params;
                if isfield(cov.signals_freqs, 'cov1d')
                    cov.signals_freqs.cov1d_freqs = cov.signals_freqs.cov1d;
                    cov.signals_freqs = rmfield(cov.signals_freqs, 'cov1d');
                end
                if isfield(cov.signals_freqs, 'cov2d')
                    cov.signals_freqs.cov2d_freqs = cov.signals_freqs.cov2d;
                    cov.signals_freqs = rmfield(cov.signals_freqs, 'cov2d');
                end

                cov.signals_freqs.Bg = Bgf;

            end
            cov.num_samples = size(this.P.pts,2);
            cov.hash_samples = DataHash(this.P.pts);
            % for resuming coverage computation;
            % only checks pts even if we have signals...
            this.CoverRes = cov;


            %% Aux function
            function  this_cov = compute_coverage(proj)
                % computes coverage for one projection

                idim = numel(proj);
                proj_st = proj{1};
                for j = 2:idim
                    proj_st = [proj_st '__x__' proj{j}];
                end
                num_grid = 1;
                for ip =  1:numel(proj)
                    num_grid = num_grid*numel(opts.params.(proj{ip}).breach_domain.enum);
                end
                vals = Bg.GetParam(proj);
                vals = unique(vals','rows');
                this_cov = 100*size(vals,1)/num_grid;
                cov.params.(proj_st) = this_cov;
            end

        end

        function [Bgrid, Bsub] = GridFilter(this, params, delta, Bgrid)
            %  Bg = this.GridFilter(params, opts [, Bg])
            %
            %  opts is an option structure defined with SetCoverageOptions
            %  Bg is the set of grid points covered by this, with cell
            %  count. 
            %  This method works by counting samples into a hashtable,
            %  where keys are defined by their position in the grid (sparse
            %  encoding in the grid). 
            % 
             
            if nargin == 2
                if isstruct(params)   % first argument is a options structure already
                    delta = params;                               
                    params = fieldnames(delta.params)';
                elseif isnumeric(params)                     
                    delta=params;
                    params = this.GetParamList();                   
                end
            end
           
            if ~exist('params', 'var')||isequal(params, 'all')
                params = this.GetParamList();                
            end

            params = setdiff(params,{'traj__ref__'}, 'stable'); 

            if  ~exist('delta','var')||isempty(delta)
                delta = this.SetCoverageOptions();
            end
            
            if isnumeric(delta)
                if isscalar(delta)
                    delta= ones(1, numel(params))*delta;
                end
                for id= 1:numel(params)
                        opts.params.(params{id}).grid = delta(id);                        
                end                
                opts = this.SetCoverageOptions(opts);
            elseif iscell(delta)
                for id= 1:numel(params)
                        opts.params.(params{id}).grid = delta{id};                        
                end
                 opts = this.SetCoverageOptions(opts);
            else  % supports legacy delta, scalar or array to specify grid 
                opts = this.SetCoverageOptions(delta);
            end            
            
            if isempty(params)
                params = fieldnames(opts.params)';
            end

            % gets grid specs for params
            grids = cell(1,numel(params));
            sizes = cell(1,numel(params));
            for ip = 1:numel(params)
                param_ip = params{ip};
                grids{ip} = opts.params.(param_ip).grid;
                sizes{ip} =  opts.params.(param_ip).eps;
                if isfield(opts.params.(param_ip), 'breach_domain')
                    grid_domain(ip) = opts.params.(param_ip).breach_domain;
                else
                    opts.params.(param_ip).breach_domain = BreachDomain('double', [opts.params.(param_ip).grid(1) opts.params.(param_ip).grid(end)]);
                    grid_domain(ip) = opts.params.(param_ip).breach_domain;
                end
            end
            % TODO: double check this for iterative computation
            if ~isempty(this.CoverRes)
                last_idx_in  = this.CoverRes.num_samples;
            else
                last_idx_in = 0;
            end

            pts_idx =last_idx_in+1:size(this.P.pts, 2);
            pts = this.GetParam(params,last_idx_in+1:size(this.P.pts, 2));
            bg_params= [params {'count' 'idx'}];
            nb_params=  numel(params);
            traj_ref = this.GetParam('traj__ref__', pts_idx); % returns [] if no traj ref
            
            % TODO that too
            if ~exist('Bg','Var')||isempty(Bgrid)
                Bgrid = BreachSet(bg_params);
                Bgrid.SetDomain(bg_params, [grid_domain BreachDomain('int'), BreachDomain('int')]);
            end

            % From now on we have Bgrid, a breach set with parameters
            % params, subset of parameters of  *this* 

            for i_pts = 1:size(pts,2)
                keep = 1; % or not
                pt = pts(:,i_pts);
                center_pt =0* pt;
                epsi = 0*pt;
                for i_dim=1:numel(params)
                    % if outside range, throw away
                    if pt(i_dim)<grid_domain(i_dim).domain(1)||...
                       pt(i_dim)>grid_domain(i_dim).domain(2)
                       keep = 0; % out
                       break;
                    end 

                    cell_centers = grid_domain(i_dim).enum;
                    [~ , imin] = min(abs(cell_centers-pt(i_dim))); % finds closest element
                    center_pt(i_dim,1) = cell_centers(imin);
                    epsi(i_dim) = sizes{i_dim}(imin);
                end
                if keep==0
                    continue;
                end

                % center_pt is now a grid pts, we add it to the map
                hash = DataHash(center_pt);
                if Bgrid.DeltaGridMapObj.isKey(hash)    
                    grid_pt = Bgrid.DeltaGridMapObj(hash); % grid pt is in already
                    if isempty(traj_ref) % we count samples
                        grid_pt.count = grid_pt.count+1;
                        grid_pt.idx(end+1) = i_pts;
                        grid_pt.dist2grid(end+1) = sqrt(sum((pt-center_pt).^2));  % distance from point to center, see if L2 dist should be sth else ..?
                    else  % we count traces
                        i_traj = traj_ref(i_pts);
                        grid_pt.idx = union(grid_pt.idx, i_traj,'stable'); % add unique
                        grid_pt.count = numel(grid_pt.idx);
                        grid_pt.dist2grid(end+1) = sqrt(sum((pt-center_pt).^2));  % discrepancy between idx and dist2grid, is this bad ?
                    end
                else
                    grid_pt.pt = center_pt;
                    grid_pt.epsi = epsi;
                    grid_pt.dist2grid = sqrt(sum((pt-center_pt).^2));  % distance from point to center, see if L2 dist should be sth else ..? 
                    if isempty(traj_ref)
                        grid_pt.count = 1;
                        grid_pt.idx = i_pts;
                    else % we count traces
                        i_traj = traj_ref(i_pts);
                        grid_pt.idx = i_traj;
                        grid_pt.count = 1;
                    end
                end
                Bgrid.DeltaGridMapObj(hash) = grid_pt;
         
            end
            grid_pts = Bgrid.DeltaGridMapObj.values;
            num_pts = numel(grid_pts);
            Bgrid.P.dim = 1:nb_params;
            Bgrid.P.epsi = zeros(nb_params, num_pts);
            
            % Legacy SetParam...
            if isempty(traj_ref) % we count samples
                for i_pt= 1:num_pts
                    [~, idx_min_dist] = min(grid_pts{i_pt}.dist2grid);

                    Bgrid.P.pts(:,i_pt) = [grid_pts{i_pt}.pt;
                        grid_pts{i_pt}.count;
                        grid_pts{i_pt}.idx(idx_min_dist)];

                    Bgrid.P.epsi(:,i_pt) = grid_pts{i_pt}.epsi;
                end
            else   % for traces, see how we can still exploit min_dist to centers...
                for i_pt= 1:num_pts

                    Bgrid.P.pts(:,i_pt) = [grid_pts{i_pt}.pt;
                        grid_pts{i_pt}.count;
                        grid_pts{i_pt}.idx(1)];

                    Bgrid.P.epsi(:,i_pt) = grid_pts{i_pt}.epsi;
                end

                    
                    end
            
            Bgrid.P = Preset_traj_ref(Bgrid.P);

            if nargout>=2
                idx_sub = Bgrid.GetParam('idx');
                Bsub=  this.ExtractSubset(idx_sub);
            end


            function delta = checks_delta(n, delta)
                % first checks if delta is all integers >=1

                % checks if cells
                if iscell(delta)
                    if numel(delta)~=n
                        if numel(delta) == 1
                            delta= repmat(delta, 1,n);
                        else
                            error('argument is not a valid grid specification, wrong size');
                        end
                    end
                    % nothing more.
                else
                    % not a cell
                    if ~isnumeric(delta) % what then ?
                        error('argument is not a valid grid specification, should be a scalar, and array or a cell');
                    end

                    if all(delta==floor(delta))&&all(delta>=1) % numbers
                        if isscalar(delta)
                            delta= repmat({delta}, 1,n);
                        elseif size(delta,1)==1||size(delta, 2) == 1
                            if numel(delta) ~= n
                                error('argument is not a valid grid specification, wrong size');
                            else
                                delta = num2cell(delta);
                            end
                        end

                    elseif size(delta,1)==1||size(delta, 2) == 1
                        delta = repmat({delta}, 1, n);
                    elseif size(delta,1)== n
                        for i = 1:n
                            delta_out{i} = delta(i, :);
                        end
                        delta= delta_out;
                    elseif size(delta,2)== n
                        for i = 1:n
                            delta_out{i} = delta(i,:);
                        end
                        delta = delta_out;
                    end
                end
            end

        end


        function [Bg,Bf,Ba] = GridFilterSignals(this, signals,  varargin)
       % [Bg, Bf, Ba] = GridFilterSignals(this, signals,  varargin)
       % 
       % - signals is the list of signals to project on, including (or not),
       % time. 
       % - varargin is the same as SetCoverageOpts
       % 
       %  Returns filtered data, grid data and all.
            
            switch nargin
                case 1  % no argument, include all signals with default grid spec (10 bins each dim)
                    sigs = this.GetSignalList();
                    params = ['time' sigs]; % parameters for grid filtering call           
                case 2  % signals is specified or a struct is provided
                    if isstruct(signals)
                        params = ['time' this.GetSignalList()];
                        opts.params = signals.signals; % assumes a valid options structure...  
                        varargin = {opts};
                    else                        
                        params =signals; % time has to be explicitly requested here                    
                    end                
                otherwise
                    if ischar(signals)  %  allow for a pair of options arguments to be started - probably confusing
                                                %  yup, I'm confused.
                        params = ['time' this.GetSignalList()];
                        varargin = [{signals} varargin]; 
                    else 
                        params =signals;   
                        if isstruct(varargin{1})
                            opts_sigs = varargin{1};
                            opts.params = opts_sigs.signals; % assumes a valid options structure...
                            varargin{1} = opts;
                        end
                    end                    
            end
                       
            %% fill in signal values
            all_sigs = this.GetSignalList();            
            all_params = ['time' all_sigs]; % all params for grid filtering
            Ba = BreachSet(all_params); 
            
            sig_values = this.GetSignalValues(all_sigs);
            if ~iscell(sig_values)
                sig_values={sig_values};
            end
            
            for itraj = 1:numel(sig_values)
                this_time= this.P.traj{itraj}.time;
                val =[this_time ; sig_values{itraj}];
                % cleaning nan in trace
                idx_not_nan = find(sum(isnan(val), 1)==0);
                val = val(:,idx_not_nan);

                if itraj==1
                    Ba.SetParam(all_params, val);
                else
                    Ba.AddParam(all_params, val);
                end
            end

            %% Applies grid filter                     
           [Bg,Bf] = Ba.GridFilter(params, varargin{:}); 
            

        end

 

        function idx_in_this = get_cell_content_pts(this, param, value)
            if isempty(this.CoverRes)
                error('Compute coverage first with SetCoverageOpts and GetCoverage.')
            end
            Bg  = this.CoverRes.params.Bg;
            domains = Bg.GetDomain(param);
            for ip = 1:numel(domains)
                value(ip) = domains(ip).checkin(value(ip));
            end
            all_pts = Bg.GetParam(param);
            idx_val_in_Bg = find(all(bsxfun(@eq, all_pts, value), 1)); % find value in Bg.P.pts
            idx_in_this = [];
            for ip = idx_val_in_Bg
                pt = Bg.P.pts(1:end-2, ip);
                hash_pt = DataHash(pt);
                pt_info = Bg.DeltaGridMapObj(hash_pt);
                idx_in_this = [idx_in_this pt_info.idx];
            end
        end

        function idx_in_this = get_cell_content_signals(this, param, value)
            if isempty(this.CoverRes)
                error('Compute coverage first with SetCoverageOpts and GetCoverage.')
            end
            Bg  = this.CoverRes.signals.Bg;
            domains = Bg.GetDomain(param);
            for ip = 1:numel(domains)
                value(ip) = domains(ip).checkin(value(ip));
            end
            all_pts = Bg.GetParam(param);
            idx_val_in_Bg = find(all(bsxfun(@eq, all_pts, value), 1)); % find value in Bg.P.pts
            idx_in_this = [];
            for ip = idx_val_in_Bg
                pt = Bg.P.pts(1:end-2, ip);
                hash_pt = DataHash(pt);
                pt_info = Bg.DeltaGridMapObj(hash_pt);
                idx_in_this = [idx_in_this pt_info.idx];
            end
        end

        function idx_in_this = get_cell_content_freqs(this, param, value)
            if isempty(this.CoverRes)
                error('Compute coverage first with SetCoverageOpts and GetCoverage.')
            end
            Bg  = this.CoverRes.signals_freqs.Bg;
            domains = Bg.GetDomain(param);
            for ip = 1:numel(domains)
                value(ip) = domains(ip).checkin(value(ip));
            end
            all_pts = Bg.GetParam(param);
            idx_val_in_Bg = find(all(bsxfun(@eq, all_pts, value), 1)); % find value in Bg.P.pts
            idx_in_this = [];
            for ip = idx_val_in_Bg
                pt = Bg.P.pts(1:end-2, ip);
                hash_pt = DataHash(pt);
                pt_info = Bg.DeltaGridMapObj(hash_pt);
                idx_in_this = [idx_in_this pt_info.idx];
            end
        end

        function idx_traces = get_traces_in_cell(this, param,value)
            if isempty(this.CoverRes)
                error('Compute coverage first with SetCoverageOpts and GetCoverage.')
            end
            if size(value, 2)>1
                value=value';
            end

            Bg  = this.CoverRes.signals.Bg;
            domains = Bg.GetDomain(param);
            for ip = 1:numel(domains)
                value(ip) = domains(ip).checkin(value(ip));
            end
            all_pts = Bg.GetParam(param);
            idx_val_in_Bg = find(all(bsxfun(@eq, all_pts, value), 1)); % find value in Bg.P.pts
            idx_traces = [];
            for ip = idx_val_in_Bg
                pt = Bg.P.pts(1:end-2, ip);
                hash_pt = DataHash(pt);
                pt_info = Bg.DeltaGridMapObj(hash_pt);
                idx_traces = [idx_traces pt_info.idx];
            end
        end

        function idx_traces = get_min_cover_traces(this)

            if isempty(this.CoverRes)
                error('Compute coverage first with SetCoverageOpts and GetCoverage.')
            end

            Bg  = this.CoverRes.signals.Bg;
            all_pts_info = Bg.DeltaGridMapObj.values;

            for ip = 1:numel(all_pts_info)
                num_traces(ip) = all_pts_info{ip}.count;
            end
            % sort by number of traces in each cell, first ones only contain
            % one trace
            [~, ia] = sort(num_traces);
            idx_traces = [];   % current traces for coverage
            for ii = ia % loop over cells in this order
                traces_ia = all_pts_info{ii}.idx;  % traces contained in this cell
                if isempty(intersect(traces_ia, idx_traces)) % is there no trace in current traces that cover this cell ?
                    idx_traces = [idx_traces  traces_ia(1)];    % we pick the first one. This is where we could optimize, pick the most covering trace
                end
            end

        end

        function idx_traces = get_min_cover_traces_freqs(this)

            if isempty(this.CoverRes)
                error('Compute coverage first with SetCoverageOpts and GetCoverage.')
            end

            Bg  = this.CoverRes.signals_freqs.Bg;
            all_pts_info = Bg.DeltaGridMapObj.values;

            for ip = 1:numel(all_pts_info)
                num_traces(ip) = all_pts_info{ip}.count;
            end
            % sort by number of traces in each cell, first ones only contain
            % one trace
            [~, ia] = sort(num_traces);
            idx_traces = [];   % current traces for coverage
            for ii = ia % loop over cells in this order
                traces_ia = all_pts_info{ii}.idx;  % traces contained in this cell
                if isempty(intersect(traces_ia, idx_traces)) % is there no trace in current traces that cover this cell ?
                    idx_traces = [idx_traces  traces_ia(1)];    % we pick the first one. This is where we could optimize, pick the most covering trace
                end
            end

        end

        function is_covered = check_covered_pts(this, params, pts)
        % checks if pts are covered by the current grid set - grid set
        % needs to be computed already using 
        % FIXME (?): assumes this is a grid set, with enum domains        


        is_covered = zeros(1, size(pts, 2));
        [iparams, ifound] = FindParam(this.P, params);

        if all(ifound)
            pts_in_this = this.GetParam(params);
            for ipt = 1:size(pts,2)
                pt  = pts(:,ipt);
                is_in = 1; % is inside the domain or not ?

                for idim =1:numel(iparams)
                    this_idx_param = iparams(idim);

                    if this.Domains(this_idx_param).is_out(pt(idim))
                        is_in = 0; % not inside the domain
                        break                        
                    else
                        pt(idim) = this.Domains(this_idx_param).checkin(pt(idim));
                    end
                end

                if is_in
                    % pt is now a grid pt if grid there is
                    if ~isempty(this.DeltaGridMapObj) % checks in map if not empty
                        hash_pt = DataHash(pt);
                        if this.DeltaGridMapObj.isKey(hash_pt);
                            is_covered(ipt)=1;
                        else
                            idx_pt_in_this = find(all(bsxfun(@eq, pts_in_this, pt)), 1);
                            if ~isempty(idx_pt_in_this)
                                is_covered(ipt)=1;
                            end
                        end
                    end
                end
            end
        else
            not_found = find(~ifound);
            warning('Parameter %s not found.', params(not_found));
        end


        end

        function h = plot_cover_stats(this)
            % Plots coverage results in all projections as percentage bars
            if ~isempty(this.CoverRes)
                cov = this.CoverRes;
                if isfield(cov,'params')
                    all_stats = fieldnames(cov.params)';
                    labels = {};
                    for istat = 1:numel(all_stats)
                        stat = all_stats{istat};
                        if strcmp(stat,'Bg')||~isempty(regexp(stat, 'cov[0-9]+d'))
                            continue;
                        end

                        params_stat = strsplit(stat,'__x__');
                        num_dim  = numel(params_stat);
                        % projection of size num_dim, if not there
                        % already, creates it and adds projection label
                        if num_dim>size(labels, 2)
                            labels{1,num_dim} ={};
                            max_dim = num_dim;
                        end
                        labels{num_dim}{end+1} = stat;
                    end

                    % loops over all dimensions:
                    h = [];
                    for dim = 1:max_dim
                        covxd = ['cov' num2str(dim) 'd'];
                        if isfield(cov.params,covxd)
                            h(end+1) = figure;
                            % first thing, get all fields of struct, barh the values
                            labels_xd = labels{dim};
                            values_xd = cellfun(@(name)(cov.params.(name)), labels_xd);
                            b = barh(values_xd);
                            yticklabels(labels_xd);
                            set(gca, 'TickLabelInterpreter', 'None', 'XLim', [0 100])
                            title([num2str(dim) 'd param coverage: ' num2str(cov.params.(covxd)) ' %']);
                            grid on;
                        end
                    end
                end

                if isfield(cov,'signals')
                    all_stats = fieldnames(cov.signals)';
                    labels = {};
                    time_labels = {};
                    for istat = 1:numel(all_stats)
                        stat = all_stats{istat};
                        % is it not a projection
                        if any(strcmp(stat,{'Bg', 'Bs_all'}))||~isempty(regexp(stat, 'cov[0-9]+d'))
                            continue;
                        end

                        signals_stat = strsplit(stat,'__x__');
                        % is it time ?

                        if ismember('time', signals_stat)
                            num_dim  = numel(signals_stat)-1;
                            % projection of size num_dim, if not there
                            % already, creates it and adds projection label
                            if num_dim>size(time_labels, 2)
                                time_labels{1,num_dim} ={};
                                max_dim_time = num_dim;
                            end
                            time_labels{num_dim}{end+1} = stat;
                        else % it is not time
                            num_dim  = numel(signals_stat);
                            % projection of size num_dim, if not there
                            % already, creates it and adds projection label
                            if num_dim>size(labels, 2)
                                labels{1,num_dim} ={};
                                max_dim = num_dim;
                            end
                            labels{num_dim}{end+1} = stat;
                        end
                    end
                    % loops over all dimensions:

                    h = [];
                    for dim = 1:max_dim_time
                        covxd = ['cov' num2str(dim) 'd_time_sig'];
                        if isfield(cov.signals,covxd)
                            h(end+1) = figure;
                            % first thing, get all fields of struct, barh the values
                            labels_xd = labels{dim};
                            values_xd = cellfun(@(name)(cov.signals.(name)), labels_xd);
                            b = barh(values_xd);
                            yticklabels(labels_xd);
                            set(gca, 'TickLabelInterpreter', 'None', 'XLim', [0 100])
                            title([num2str(dim) 'd time signals coverage: ' num2str(cov.signals.(covxd)) ' %']);
                            grid on;
                        end
                    end

                    for dim = 1:max_dim
                        covxd = ['cov' num2str(dim) 'd_sig'];
                        if isfield(cov.signals,covxd)
                            h(end+1) = figure;
                            % first thing, get all fields of struct, barh the values
                            labels_xd = labels{dim};
                            values_xd = cellfun(@(name)(cov.signals.(name)), labels_xd);
                            b = barh(values_xd);
                            yticklabels(labels_xd);
                            set(gca, 'TickLabelInterpreter', 'None', 'XLim', [0 100])
                            title([num2str(dim) 'd signals coverage: ' num2str(cov.signals.(covxd)) ' %']);
                            grid on;
                        end
                    end
                end

                if isfield(cov,'signals_freqs')
                    all_stats = fieldnames(cov.signals_freqs)';
                    labels = {};
                    for istat = 1:numel(all_stats)
                        stat = all_stats{istat};
                        if any(strcmp(stat,{'Bg', 'Bs_all'}))||~isempty(regexp(stat, 'cov[0-9]+d'))
                            continue;
                        end

                        params_stat = strsplit(stat,'__x__');
                        num_dim  = numel(params_stat);
                        % projection of size num_dim, if not there
                        % already, creates it and adds projection label
                        if num_dim>size(labels, 2)
                            labels{1,num_dim} ={};
                            max_dim = num_dim;
                        end
                        labels{num_dim}{end+1} = stat;
                    end

                    % loops over all dimensions:
                    h = [];
                    for dim = 1:max_dim
                        covxd = ['cov' num2str(dim) 'd'];
                        if isfield(cov.signals_freqs,covxd)
                            h(end+1) = figure;
                            % first thing, get all fields of struct, barh the values
                            labels_xd = labels{dim};
                            values_xd = cellfun(@(name)(cov.signals_freqs.(name)), labels_xd);
                            b = barh(values_xd);
                            yticklabels(labels_xd);
                            set(gca, 'TickLabelInterpreter', 'None', 'XLim', [0 100])
                            title([num2str(dim) 'd frequencies coverage: ' num2str(cov.signals_freqs.(covxd)) ' %']);
                            grid on;
                        end
                    end
                end

            end
        end

        function h = plot_cover_proj(this,proj, cov_type)
            
            if ~isempty(this.CoverRes)
                cov   = this.CoverRes;
                opts = this.CoverOpts;

                % find projections
                if nargin<3
                    cov_type =[];
                end

                % look into parameters
                if isfield(cov,'params')
                    if isfield(cov.params,proj)
                        Bg= cov.params.Bg;
                        cov_type = 'params';
                    end
                end

                % look into signals, if not found yet
                if isfield(cov,'signals')&&isempty(cov_type)
                    if isfield(cov.signals,proj)
                        Bg= cov.signals.Bg;
                        cov_type = 'signals';
                    end
                end

                % look into freqs, if not found yet
                if isfield(cov,'signals_freqs')&&isempty(cov_type)
                    if isfield(cov.signals_freqs,proj)
                        Bg= cov.signals_freqs.Bg;
                        cov_type = 'signals_freqs';
                    end
                end
                if isempty(cov_type)
                    error('plot_cover_proj:proj_not_found',['Projection ' proj ' not found.']);
                end

                h = gcf;
                params= strsplit(proj,'__x__'); % assumes p1__x__p2__x__...
                if numel(params)>2
                    warning('more than two parameters, not supported');
                    return;
                end

                values= Bg.GetParam(params);

                [vu, ~, iv] = unique(values','row');
                vu = vu'; 
                iv=iv';
                switch cov_type
                    case 'params'
                        for idx = 1:size(vu,2) % loop over cells that share the same projections
                            idx_in_set  = this.get_cell_content_pts(params, vu(:,idx));
                            cu(idx) = numel(unique(idx_in_set));
                        end
                    case 'signals'
                        for idx = 1:size(vu,2) % loop over cells that share the same projections
                            idx_in_set  = this.get_cell_content_signals(params, vu(:,idx));
                            cu(idx) = numel(unique(idx_in_set));
                        end
                    case 'signals_freqs'
                        for idx = 1:size(vu,2) % loop over cells that share the same projections
                            idx_in_set  = this.get_cell_content_freqs(params, vu(:,idx));
                            cu(idx) = numel(unique(idx_in_set));
                        end
                end

                switch(numel(params))
                    case 1
                        x_all = opts.(cov_type).(params{1}).breach_domain.enum;
                        widthx = opts.(cov_type).(params{1}).eps; % TODO replace or update bar function call to use width of bars
                        [~,~,idx_vu] = intersect(vu, x_all); % non empty bins

                        y = zeros(1,numel(x_all));
                        y(idx_vu) = cu;

                        bar(x_all, y);
                        hold on;
                        grid on;

                    case 2

                        x_all = opts.(cov_type).(params{1}).breach_domain.enum;
                        widthx = 2*opts.(cov_type).(params{1}).eps;
                        y_all = opts.(cov_type).(params{2}).breach_domain.enum;
                        widthy = 2*opts.(cov_type).(params{2}).eps;

                        x = vu(1,:);
                        y = vu(2,:);
                        plot3(x,y,cu, 's');
                        hold on;
                        view(0,90);
                        cu_matrix = zeros(numel(y_all),numel(x_all));
                        for i_x =1:numel(x_all)
                            idx = find(x== x_all(i_x));
                            [~, ~, idx_y] = intersect(y(idx), y_all);
                            cu_matrix(idx_y,i_x) = cu(idx);
                        end
                        [X, Y] = meshgrid(x_all,y_all);

                        % colormap, copper reversed
                        cmap= flipud(colormap('gray'));
                        cmap = [1 1 1;
                            repmat(cmap(32,:),32,1);
                            cmap(33:220,:);
                            repmat(cmap(220,:),46,1)];
                        colormap(cmap);

                        grid on;
                        %scatterbar3(X,Y,cu_matrix,widthx, widthy); % TODO: option for 3d plot ?
                        flatbar3(X,Y,cu_matrix,widthx, widthy);
                        caxis([0 max(max(cu),1e-3)]);
                        c = colorbar;
                        ticks = get(c,'Ticks');
                        ticks = unique([0 floor(ticks(ticks>0))], 'stable');
                        set(c, 'Ticks',ticks);

                        xlabel(params{1},'Interpreter','none');
                        ylabel(params{2},'Interpreter','none');
                        zlabel('#samples');



                    otherwise
                        return;
                end
                title([proj ' ' cov_type ' coverage: ' num2str(cov.(cov_type).(proj)) ' %'],'interpreter', 'none');
            end

        end
               
        function Bg_proj = SimpleProject(this,params)
        % Projects on provided dimensions
            Bg_proj = BreachSet(params);
            values= this.GetParam(params);
            [vu, iv] = unique(values','row');
            Bg_proj.SetParam(params,vu');
            Bg_proj.P.dim = 1:numel(params);
            Bg_proj.P.epsi = zeros(numel(params), size(vu, 1));

            idx_params = FindParam(this.P, params);          
            for ip = 1:numel(idx_params)
                idx_this_param = idx_params(ip);
                idx_epsi= find(this.P.dim==idx_this_param,1);
                if ~isempty(idx_epsi)
                    Bg_proj.P.epsi(idx_this_param,:) = this.P.epsi(idx_epsi,iv);
                end
            end
            % we might want to remove 0 epsi from dim, if you see what I
            % mean (surely not if you're not me)
            
            
        end
       
        function r_list = plot_cover_grid(this, proj, varargin)

            Bproj = this.SimpleProject(proj);            

            P = Bproj.P;
            if ~isnumeric(proj)
                proj = FindParam(P,proj);
            end
            col = 'b';
            alph = .1;
            proj = proj(proj>0);
            proj = proj(proj<=size(P.pts,1));
            if isempty(proj)
                proj = 1:min(3,numel(P.dim));
                proj = P.dim(proj);
            elseif(numel(proj)>3)
                proj = proj(1:3);
            end

            if(~exist('ipts','var')||isempty(ipts))
                ipts = 1:size(P.pts,2);
            end

            if ~exist('col','var')
                col = 'b';
            end

            if ~exist('alph','var')
                alph = .1;
            end

            hold on;
            r_list = [];
            switch(numel(proj))

                case 1
                    if isfield(P,'ParamList')
                        xlabel(P.ParamList{proj(1)},'Interpreter','none');
                    else
                        xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
                    end

                    DX(2) = 0;
                    idx = find(P.dim == proj);

                    for kk=ipts % for each parameter sample
                        if ~isempty(idx)
                            DX(1) = P.epsi(idx,kk);
                        else
                            DX(1) = 0;
                        end
                        X = [P.pts(proj,kk)-DX(1),-DX(2)];
                        r=  rect(X, 2*DX, col, alph);                       
                        set(r, varargin{:});

                    end
                    %set(gca, 'YLim', [-1 1], 'YtickLabel', {},'Interpreter','none');

                case 2
                    if isfield(P,'ParamList')
                        xlabel(P.ParamList{proj(1)},'Interpreter','none');
                        ylabel(P.ParamList{proj(2)},'Interpreter','none');
                    else
                        xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
                        ylabel(['x_' num2str(proj(2))],'Interpreter','tex');
                    end

                    DX = zeros(1,2);
                    for kk=ipts
                        for jj = 1:2 % if proj(jj) in P.dim, set the width
                            idx = find(P.dim == proj(jj));
                            if isempty(idx)
                                DX(jj) = 0;
                            else
                                DX(jj) = P.epsi(idx,kk);
                            end
                        end
                        X = P.pts(proj,kk)'-DX;
                        r = rect(X, 2*DX, col, alph);
                        set(r, varargin{:});
                    end

                case 3
                    if isfield(P,'ParamList')
                        xlabel(P.ParamList{proj(1)},'Interpreter','none');
                        ylabel(P.ParamList{proj(2)},'Interpreter','none');
                        zlabel(P.ParamList{proj(3)},'Interpreter','none');
                    else
                        xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
                        ylabel(['x_' num2str(proj(2))],'Interpreter','tex');
                        zlabel(['x_' num2str(proj(3))],'Interpreter','tex');
                    end

                    DX = zeros(1,3);
                    for kk=ipts
                        for jj = 1:3 % if proj(jj) in P.dim, set the width
                            idx = find(P.dim == proj(jj));
                            if isempty(idx)
                                DX(jj) = 0;
                            else
                                DX(jj) = P.epsi(idx,kk);
                            end
                        end

                        X = P.pts(proj,kk)'-DX;
                        r= voxel(X, 2*DX, col, alph);
                        set(r, varargin{:});
                    end
            end
            r_list(end+1) =r;

            grid on;
            hold off;
            
        end

        function Bp  = sig2params(this)

            sigs = this.GetSignalList();
            params = ['traj__ref__' 'time' sigs];
            Bp = BreachSet(params);

            sig_values = this.GetSignalValues(sigs);
            if ~iscell(sig_values)
                sig_values={sig_values};
            end

            for itraj = 1:numel(sig_values)
                val = [ones(1,size(sig_values{itraj},2))*itraj; this.P.traj{itraj}.time; sig_values{itraj}]; % keep itraj for traj_ref__
                if itraj==1
                    Bp.SetParam(params, val);
                else
                    Bp.AddParam(params, val);
                end
            end

            X = Bp.GetParam(sigs);
            minX = min(X,[],2);
            maxX = max(X,[],2);
            tolAmp = (maxX-minX)/200;
            Bp.SetParamRanges(sigs, [(minX - tolAmp), (maxX+tolAmp)]);

        end


        %% Requirements - sort of deprecated per BreachRequirement class
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
            options = varargin2struct_breach(options, varargin{:});

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
            options = struct('FolderName', []);
            options = varargin2struct_breach(options, varargin{:});

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
                        traces(it).specs.params(ip).names = fieldnames(params)';
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
            options = varargin2struct_breach(options, varargin{:});

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

        %% Signatures
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

            f = fieldnames(sigp)';
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

                if idx==inf
                    warning('BreachSet:GetSignalSignature:not_found', 'Signal or alias %s not found.', sig);
                end

                if idx==is % first time we see this guy, take it as rep
                    sigs.signals_reps{end+1} = sig;
                    idx = numel(sigs.signals_reps);
                    sigs.signals_map_idx(is) = idx;
                else % idx is smaller than is so signals_map_idx(idx) is defined and correct, not necessarilly equal to idx ...
                    sigs.signals_map_idx(is) = sigs.signals_map_idx(idx);
                end


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

        function varargout = disp_signals(this)
            max_length = 200;
            sig_list = this.P.ParamList(1:this.P.DimX);
            if isfield(this.P, 'traj')
                num_traces = numel(this.P.traj);
            else
                num_traces = 0;
            end
         
            switch num_traces
                case 0
                    str_num = '0 trace.\n';
                case 1
                    str_num = '1 trace.\n';
                otherwise
                    str_num= [num2str(num_traces) ' traces\n'];
            end                      

            switch numel(sig_list)
                case 0
                    str = '0 signals.\n';
                case 1
                    str = ['1 signal: ' sig_list{1} '.\n' str_num];
                otherwise
                    str_sig_list = list_manip.to_string(sig_list, ', ');
                    if numel(str_sig_list) > max_length
                        str_sig_list = [str_sig_list(1:max_length) ', ...'];
                    end
                    str = [num2str(numel(sig_list)) ' signals: ' str_sig_list '.\n'  str_num];
            end
            if nargout == 0
                varargout = {};
                fprintf(str);
            else
                varargout{1} = str;
            end
        end

        function varargout = disp_params(this)
            max_length = 200;
            params_list = this.P.ParamList(this.P.DimX+1:end);
            num_samples = this.GetNbParamVectors();

            switch num_samples
                case 0
                    str_num = '0 samples.\n';
                case 1
                    str_num = '1 sample.\n';
                otherwise
                    str_num= [num2str(num_samples) ' samples\n'];
            end

            switch numel(params_list)
                case 0
                    str = '0 Parameters.\n';                    
                case 1
                    str = ['One Parameter: ' params_list{1} '.\n' str_num];
                otherwise
                    str_params_list = list_manip.to_string(params_list, ', ');
                    if numel(str_params_list) > max_length
                        str_params_list = [str_params_list(1:max_length) ', ...'];
                    end
                    str = [num2str(numel(params_list)) ' Parameters: ' str_params_list '.\n' str_num];
            end



            if nargout == 0
                varargout = {};
                fprintf(str);
            else
                varargout{1} = str;
            end
        end

        function varargout = disp_short_description(this)
            
            param_desc = this.disp_params();

            str = ['----\n' ...
                      param_desc ...
                      '----\n'];

            if this.P.DimX>0
                sig_desc = this.disp_signals();
                str= [str ...
                    sig_desc ...
                    '----\n'];
            end

            if nargout == 0
                varargout = {};
                fprintf(str);
            else
                varargout{1} = str;
            end

        end


        function varargout = disp(this)
        
            str = 'BreachSet object.\n';
            str = [str this.disp_short_description()]; 
            
            if nargout == 0
                varargout = {};
                fprintf(str);
            else
                varargout{1} = str;
            end
            
        end
        


        %% Stuff

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

        function varargout = PrintSignals(this)
            st = sprintf('---- SIGNALS ----\n');
            for isig = 1:this.P.DimX
                st = sprintf([st '%s %s\n'], this.P.ParamList{isig}, this.get_signal_attributes_string(this.P.ParamList{isig}));
            end
            st = sprintf([st '\n']);
            st = [st this.PrintAliases()];

            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end

        end

        function varargout = PrintAliases(this)
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
            end
            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
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

        function varargout = PrintParams(this,params, header)
            st = '';
            nb_pts= this.GetNbParamVectors();

            if ~exist('params','var')
                params = this.P.DimX+1:numel(this.P.ParamList);
            elseif ischar(params)||iscell(params)
                params = FindParam(this.P, params);
            end

            if ~exist('header', 'var')
                header = '-- PARAMETERS --';
            end

            if (nb_pts<=1)
                if ~isempty(header)
                    st = [header sprintf('\n')];
                end
                for ip = params
                    st = [st sprintf('%s=%g       %s\n',this.P.ParamList{ip},this.P.pts(ip,1), this.Domains(ip).short_disp(1))];
                end

            else
                if ~isempty(header)
                    st = [header sprintf(' (%d samples):\n',nb_pts)];
                end
                for ip = params
                    st = [st sprintf('%s     %s\n',this.P.ParamList{ip}, this.Domains(ip).short_disp(1))];
                end

            end

            st = [st sprintf('\n')];

            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end

        end

        function varargout = PrintAll(this)
            this.UpdateSignalRanges();
            st = this.PrintSignals();
            st = sprintf([st  this.PrintParams()]);

            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end
        end

    
        function s= isSignal(this,params)
            idx_s = FindParam(this.P, params);
            s = idx_s <= this.P.DimX;
        end

        % Warning handler
        function WarningResetP(this, fname)
            if this.GetNbParamVectors()>1
                warning('BreachSet:warning_will_reset_pts',['This set contains more than one parameter sample or traces - the function ' fname ' will likely erase them.'])
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
            p_req  = this.P.pts(this.P.DimP+1:end,1);

            this.P = CreateParamSet(this.P);
            try
                this.P.pts = this.Sys.p;
                if ~isempty(p_req)
                    this.P.pts = [this.P.pts; p_req];
                end
                % Add property params
                props = this.Specs.values;
                for i=1:numel(props)
                    phi = props{i};
                    params_prop = get_params(phi);
                    this.SetParamSpec(fieldnames(params_prop)', cellfun(@(c) (params_prop.(c)), fieldnames(params_prop)'),1);
                end
            end
            for ip = 1:numel(this.P.ParamList)
                this.Domains(ip).domain=[];
            end

            this.resetStatus();
        end

        function ResetEpsi(this)
            % (Legacy) Set param ranges around individual parameter samples to zero
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

            if size(this.P.pts,2)>1
                % get rid of redundant pts
                [~, iu] = unique(this.P.pts','rows');
                this.P = Sselect(this.P, iu);
            end
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
            for is = 1:numel(signals)
                sig = signals{is} ;
                if this.AliasMap.isKey(sig)
                    aliases = union(aliases, this.AliasMap(sig),'stable');
                end
            end

        end

    end


%% Old coverage code
%  The following commented code is for future potential coverage measure,
%  courtesy of Jim K. 
%
%         function coverage = ComputeLogCellOccupancyCoverage(this)
%             % Compute the log of the cell occupancy
%             % for the parameter space
% 
%             % Total number of cells in the parameter space
%             total_cells = this.TotalCellCount();
%             % Next, obtain the total number of populated cells
%             pop_cells = this.NumPoints();
%             coverage = log(pop_cells)/log(total_cells);
%         end
% 
%         function coverage = ComputeCellOccupancyCoverage(this)
%             % Compute the cell occupancy for the parameter space
% 
%             % Total number of cells in the parameter space
%             total_cells = this.TotalCellCount();
%             % Next, obtain the total number of populated cells
%             pop_cells = this.NumPoints();
%             coverage = pop_cells/total_cells;
%         end
% 
%         function coverage = ComputeEntropyCoverage(this)
%             % Compute the combinatorial entropy (also known as the
%             % coarse-grained Boltzmann entropy).
% 
%             % The entropy measure is given by the following:
%             %
%             % G(D) = p!/(n_1!\cdots n_c!)
%             %
%             % S_B(D) = log(G(D)),
%             %
%             % where p is the total number of points, and n_i is the number
%             % of points in cell i.
% 
%             % Collect the number of points in each cell and compute the
%             % denominator of G(D). At the same time, compute the total
%             % number of points p.
%             delta_cell_labels = this.DeltaGridMapObj.keys;
%             numpoints = 0;
%             GD_denominator = 1;
%             for ind = 1:length(delta_cell_labels)
%                 n_i = double(this.DeltaGridMapObj(delta_cell_labels{ind}));
%                 numpoints = numpoints + n_i;
%                 GD_denominator = GD_denominator*(factorial(n_i));
%             end
% 
%             GD_numerator = factorial(numpoints);
% 
%             GD = GD_numerator/GD_denominator;
% 
%             coverage = log(GD);
%         end



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
                cmp.addStatus(10,'The two sets have different number of parameter samples.');
                return;
            end

            % Checks system parameters
            rg_sys = this.P.DimX+1:this.P.DimP; % at this point, this is the same as other
            sys_pts_this = this.P.pts(rg_sys,:);
            sys_pts_other = other.P.pts(rg_sys,:);

            diff_sys_pts = norm(sys_pts_this-sys_pts_other);
            if (diff_sys_pts == 0)
                cmp.addStatus(0, 'The two sets have the same system parameter samples.');
            else
                cmp.addStatus(1, ['Distance between the two system parameter samples: ' num2str(diff_sys_pts)])
            end

            % Checks spec. parameters
            rg_pspec = this.P.DimP+1:size(this.P.pts,1); % at this point, this is the same as other
            if ~isempty(rg_pspec)
                spec_pts_this = this.P.pts(rg_pspec,:);
                spec_pts_other = other.P.pts(rg_pspec,:);

                diff_spec_pts = norm(spec_pts_this-spec_pts_other);
                if (diff_spec_pts == 0)
                    cmp.addStatus(0, 'The two objects have the same spec. parameter samples.');
                else
                    cmp.addStatus(1, ['Distance between the two spec. parameter samples: ' num2str(diff_spec_pts)])
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
