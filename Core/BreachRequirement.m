classdef BreachRequirement < BreachTraceSystem
    % BreachRequirement Generic requirement class - it provides helpers to get
    % data from Breach objects so different types of constraints evaluations can easily be
    % implemented.
    
    properties
        BrSet
        postprocess_signal_gens            % output generators
        req_monitors
        precond_monitors
        signals_in
        traces_vals_precond % results for individual traces & precond_monitors
        traces_vals_vac  % input vacuity, if avail
        traces_vals % results for individual traces & req_monitors
        val         % summary evaluation for all traces & req_monitors
    end
    
    properties (Access=protected)
        eval_precond_only= false  % when true, does what it says
    end
    
    methods
        
        function this = BreachRequirement(req_monitors, postprocess_signal_gens, precond_monitors)
            
            this = this@BreachTraceSystem({}, [], {'data_trace_idx_'});
            % Adds requirement monitors, at least one
            this.req_monitors = {};
            
            if nargin>=1                
                if strcmp(req_monitors,'random')
                    switch nargin 
                        case 2
                            this.AddRandReq(postprocess_signal_gens);
                            postprocess_signal_gens = {};
                        case 3
                            this.AddRandReq(postprocess_signal_gens, precond_monitors);
                            postprocess_signal_gens = {};
                            precond_monitors = {};
                        otherwise
                            this.AddRandReq();
                    end
                else                    
                    this.AddReq(req_monitors);
                end
            else
                this.AddRandReq();
            end
            
            % Add output gens
            this.postprocess_signal_gens={};
            if exist('postprocess_signal_gens', 'var')&&~isempty(postprocess_signal_gens)
                this.AddPostProcess(postprocess_signal_gens);
            end
            
            % precondition requirements
            this.precond_monitors = {};
            if  exist('precond_monitors', 'var')&&~isempty(precond_monitors)
                this.AddPreCond(precond_monitors);
            end
            
            % Reset signal map and figure out what signals are required input signals
            this.ResetSigMap();
            
        end
        
        function AddReq(this, req_monitors)
            
            [~, monitors] = get_monitors(req_monitors);
            for ifo = 1:numel(monitors)
                this.AddOutput(monitors{ifo});
                monitors{ifo}.R = this;
            end
            this.req_monitors = [this.req_monitors monitors];
            this.signals_in = this.get_signals_in();
        end
        
        function AddRandReq(this, n, Ops)
            % Adds n random stl formula with given operators                                    
            if nargin == 1
                st_phi = rand_formula(1);                
                this.AddReq(st_phi);
            else
                if nargin<2
                    n = 1;
                end
                if ~exist('Ops','var')||isempty(Ops)
                    Ops = {'and', 'ev', 'evI', 'alw'};                
                end
                st_phi = rand_formula(n,Ops);               
                this.AddReq(st_phi);
            end                        
           
        end
        
        
        function AddPostProcess(this, postprocess_signal_gens)
            
            if ~iscell(postprocess_signal_gens)
                postprocess_signal_gens = {postprocess_signal_gens};
            end
            for ifo = 1:numel(postprocess_signal_gens)
                this.AddOutput(postprocess_signal_gens{ifo});
            end
            this.postprocess_signal_gens = [this.postprocess_signal_gens postprocess_signal_gens];
            this.signals_in = this.get_signals_in();
        end
        
        function AddPreCond(this, precond_monitors)
            
            [~, monitors] = get_monitors(precond_monitors);
            for ifo = 1:numel(monitors)
                this.AddOutput(monitors{ifo});
            end
            this.precond_monitors = [this.precond_monitors monitors];
            this.signals_in = this.get_signals_in();
            
        end
        
        function this = SetSignalMap(this, varargin)
            % SetSignalMap maps signal names to signals needed for requirement
            % evaluation
            %
            % Input:  a map or a list of pairs or two cells
            
            this = SetSignalMap@BreachSet(this,varargin{:});
            
            if this.verbose >= 2
                this.PrintSigMap();
            end
            
            this.signals_in = this.get_signals_in();
            
        end
        
        function ResetSigMap(this)
            this.sigMap = containers.Map();
            this.sigMapInv = containers.Map();
            this.AliasMap =containers.Map();
            
            this.signals_in = this.get_signals_in();
        end
        
        function ResetEval(this)
            this.BrSet = [];         
            this.SetParam('data_trace_idx_', 0);
            this.P = SPurge(this.P);
            if size(this.P.pts,2)>1
               this.P.pts = unique(this.P.pts', 'rows')';
               this.P.epsi = this.P.epsi(:,1:size(this.P.pts,2));                  
            end
            this.P.selected = zeros(1, size(this.P.pts,2));
            this.P = Preset_traj_ref(this.P);
            this.traces_vals_precond = [];
            this.traces_vals_vac = [];
            this.traces_vals = [];
            this.val = [];            
        end
        
        function  [val_precond, traj_req] = evalTracePrecond(this,traj_req)
            % evalTracePrecond eval satisfaction of requirements for one trace.
            
            % eval precondition on trace
            val_precond = nan(1, numel(this.precond_monitors));
            for i_req = 1:numel(this.precond_monitors)
                req = this.precond_monitors{i_req};
                [val_precond(i_req), traj_req] = this.evalRequirement(req, traj_req);
            end
            
            if any(val_precond<0)
                traj_req.status = -1;
            end
            
        end
        
        function  [val, traj_req] = evalTrace(this,traj_req)
            % evalTrace eval satisfaction of requirements for one trace.
            
            % eval requirements on trace
            val = nan(1, numel(this.req_monitors));
            for i_req = 1:numel(this.req_monitors)
                req = this.req_monitors{i_req};
                [val(i_req), traj_req] = this.evalRequirement(req, traj_req);
            end
            
        end
        
        function traces_vals_precond = EvalPrecond(this, varargin)
            % BreachRequirement.EvalPrecond returns evaluation of the
            % requirement  preconditions
            this.eval_precond_only = true;
            
            % Collect traces from context and eval them
            [~, traces_vals_precond] = this.evalAllTraces(varargin{:});
            this.traces_vals_precond = traces_vals_precond;
            
            this.eval_precond_only = false;
        end
        
        function ResetSimulations(this)
            ResetSimulations@BreachSet(this);
            this.traces_vals =[];
            this.val = [];
            this.traces_vals_precond= [];
            this.BrSet =[];
        end
        
        function isSensitive = ...
                evalStructuralSensitivity(this, traj_index1, traj_index2, this_monitor)
            
            X1 = this.GetSignalValues(this_monitor.signals_in, traj_index1);
            X2 = this.GetSignalValues(this_monitor.signals_in, traj_index2);
            
            if all(all(X1 == X2))
                isSensitive = 0;
                return
            end
            
            idx_par_req = FindParam(this.P, this_monitor.params);
            
            % We only take p_in and time from traj1 - assume they are same
            % as traj2
            p_in = this.P.traj{traj_index1}.param(1, idx_par_req)';
            time = this.P.traj{traj_index1}.time;
            
            
            isSensitive = ...
                this_monitor.get_structural_sensitivity(time, X1, X2, p_in);
        end
        
        function [global_val, traces_vals, traces_vals_precond, traces_vals_vac] = Eval(this, varargin)
            % BreachRequirement.Eval returns evaluation of the requirement -
            % compute it for all traces available and returns min (implicit
            % conjunction)
            
            % Collect traces from context and eval them
            if isfield(this.Sys, 'use_parallel')&&this.Sys.use_parallel
                [traces_vals, traces_vals_precond, traces_vals_vac] = this.parevalAllTraces(varargin{:});
            else
                [traces_vals, traces_vals_precond, traces_vals_vac] = this.evalAllTraces(varargin{:});
            end
            % A BreachRequirement must return a single value
            global_val = min(min(traces_vals));
            global_precond_val = min(min(traces_vals_precond));
            this.val = min([global_val,-global_precond_val]);
        end
                
        function [traces_vals, traces_vals_precond, traces_vals_vac] =evalAllTraces(this,varargin)
            % BreachRequirement.evalAllTraces collect traces and apply
            % evalTrace
            new= this.getBrSet(varargin{:});
            if new
                num_eval = this.GetNbParamVectors();% numel(this.P.traj);
                traces_vals = nan(num_eval, numel(this.req_monitors));
                traces_vals_vac = nan(num_eval, numel(this.req_monitors));
                traces_vals_precond = nan(num_eval, numel(this.precond_monitors));
                
                % eval pre conditions
                if ~isempty(this.precond_monitors)
                    for it = 1:num_eval
                        for ipre = 1:numel(this.precond_monitors)
                            req = this.precond_monitors{ipre};
                            traces_vals_precond(it, ipre)  = eval_req(this,req,it);
                        end
                    end
                end
                
                % eval requirement
                execTimesForThisReq = zeros(1, numel(this.req_monitors));
                for it = 1:num_eval
                    currentTime = datestr(now, 'HH:MM:ss');
                    if num_eval > 30 && (mod(it, 10) == 0)
                        fprintf(['*** START traj ' num2str(it) '/' num2str(num_eval) ' at ' currentTime '\n']);
                    end
                    
                    if any(traces_vals_precond(it,:)<0)
                        traces_vals( it, :)  = NaN;
                    else
                        for ipre = 1:numel(this.req_monitors)
                            startTimeOfReq = tic;
                            req = this.req_monitors{ipre};
                            
                            [traces_vals(it, ipre), ...
                                traces_vals_vac(it, ipre)]  = eval_req(this,req,it);
                            execTimesForThisReq(ipre) = execTimesForThisReq(ipre) + toc(startTimeOfReq);
                        end
                    end
                    % Finished the traj
                    
                end
                
                % Display the nSlowest slowest specifications
                nSlowest = 5;
                if (numel(execTimesForThisReq) > nSlowest) && (num_eval > 30)
                    fprintf(['  Finished all trajs, ' num2str(nSlowest) ' slowest specs (out of ' num2str(numel(this.req_monitors)) '): ']);
                    [sortedExecTimes, sortedSpecIndex] = sort(execTimesForThisReq, 'descend');
                    for slowestCounter = 1:nSlowest
                        thisIndex = sortedSpecIndex(slowestCounter);
                        thisTime = sortedExecTimes(slowestCounter);
                        fprintf([this.req_monitors{thisIndex}.name ' (' num2str(thisTime) 's)']);
                        if slowestCounter < nSlowest
                            fprintf(', ');
                        end
                    end
                    fprintf('\n');
                end
                this.traces_vals_precond = traces_vals_precond;
                this.traces_vals_vac = traces_vals_vac;
                this.traces_vals = traces_vals;
            else
                traces_vals_precond = this.traces_vals_precond;
                traces_vals_vac = this.traces_vals_vac;
                traces_vals = this.traces_vals;
            end
        end
        
        function [traces_vals, traces_vals_precond, traces_vals_vac] =parevalAllTraces(this,varargin)
            % BreachRequirement.evalAllTraces collect traces and apply
            % evalTrace
            new= this.getBrSet(varargin{:});
            if new
                num_traj = numel(this.P.traj);
                traces_vals = nan(num_traj, numel(this.req_monitors));
                traces_vals_vac = nan(num_traj, numel(this.req_monitors));
                traces_vals_precond = nan(num_traj, numel(this.precond_monitors));
                
                % eval pre conditions
                num_precond_monitors = numel(this.precond_monitors);
                if ~isempty(this.precond_monitors)
                    for ipre = 1:num_precond_monitors
                        req = this.precond_monitors{ipre};
                        parfor it = 1:num_traj
                            traces_vals_precond(it, ipre)  = eval_req(this,req,it);
                        end
                    end
                end
                
                
                % eval requirement
                num_monitors = numel(this.req_monitors);
                for ipre = 1:num_monitors
                    req = this.req_monitors{ipre};                            
                    parfor it = 1:num_traj
                        if any(traces_vals_precond(it,:)<0)
                            traces_vals(it, ipre)  = NaN;
                        else
                            [traces_vals(it, ipre), traces_vals_vac(it, ipre)]  = eval_req(this,req,it);
                        end
                    end
                end
                this.traces_vals_precond = traces_vals_precond;
                this.traces_vals_vac = traces_vals_vac;
                this.traces_vals = traces_vals;
            else
                traces_vals_precond = this.traces_vals_precond;
                traces_vals_vac = this.traces_vals_vac;
                traces_vals = this.traces_vals;
            end
        end
        
        function F = PlotDiagnostics(this, idx_req_monitors, itraj)
            if nargin<2
                idx_req_monitors = 1;            
            end
            if isa(idx_req_monitors,'req_monitor')
                req_mon = {idx_req_monitors};
                idx_req_monitors = 1;
            else
                req_mon = this.req_monitors(idx_req_monitors);
            end
            
            if nargin<3
                itraj=1;
            end
            F = BreachSignalsPlot(this, {}, itraj); % empty
            for ifo =1:numel(idx_req_monitors)
                if isa(req_mon{ifo},'stl_monitor')
                    req = req_mon{ifo}; 
                    traj = this.P.traj{itraj}; 
                    t = traj.time;
                    idx_par_req = FindParam(this.P, req.params);
                    p_in = traj.param(1, idx_par_req);
                    Xin = this.GetSignalValues(req_mon{ifo}.signals_in,itraj);
                    req.init_tXp(t,Xin,p_in);
                    req.plot_diagnostics(F);
                else
                    if ~isempty(req_mon{ifo}.signals_in)
                        F.AddSignals(req_mon{ifo}.signals_in)
                    end
                    if ~isempty(req_mon{ifo}.signals)
                        F.AddSignals(req_mon{ifo}.signals)
                    end
                end
            end
            for ia = 1:numel(F.Axes)
                F.update_legend(F.Axes(ia).ax);
            end        
            
            
        end
            
        function  diag_map = Explain(this, req, trace_num) 
           % Explains returns implicant for requirement req, and trace
           % number trace_num
            
           if nargin <= 1
               req = this.req_monitors{1};
           elseif isscalar(req)
               req = this.req_monitors{req};
           end
           assert(isa(req, 'stl_monitor'));
                      
           if nargin <= 2
              trace_num = 1;  
           end
           traj = this.P.traj{trace_num}; 
           
           t = traj.time;
           idx_par_req = FindParam(this.P, req.params);
           p_in = traj.param(1, idx_par_req);
           Xin = this.GetSignalValues(req.signals_in,trace_num);           
           req.explain(t,Xin,p_in);
           diag_map = req.diag_map;
                                             
        end
                
        function h = PlotSignals(this,varargin)
            h = BreachSignalsPlot(this,varargin{:});
        end
        
        function summary = GetSummary(this, varargin)
            
            summary = GetStatement(this);
            summary.signature = this.GetSignature(varargin{:});
            summary.num_violations_per_trace =sum(this.traces_vals<0 , 2 )';
            [~,idxm ] = sort(summary.num_violations_per_trace, 2, 'descend');
            summary.idx_traces_with_most_violations = idxm;
            
        end
        
        function summary = GetStatement(this)
            
            if this.CountTraces()==1
                summary.statement = sprintf('%d trace evaluated', this.CountTraces());
            else
                summary.statement = sprintf('%d traces evaluated', this.CountTraces());
            end
            summary.num_traces_evaluated =size(this.traces_vals,1);
            summary.requirements.names = cell(1,numel(this.req_monitors));
            for ir = 1:numel(this.req_monitors)
                if isa(this.req_monitors{ir}, 'stl_monitor')
                    summary.requirements.names{ir} = this.req_monitors{ir}.formula_id;
                else
                    summary.requirements.names{ir} = class(this.req_monitors{ir});
                end
            end
            if summary.num_traces_evaluated>0
                summary.val = this.val;
                summary.requirements.rob = this.traces_vals;
                summary.requirements.rob_vac = this.traces_vals_vac;
                summary.requirements.sat = this.traces_vals >=0;
                summary.num_requirements = size(this.traces_vals,2);
                if summary.num_requirements == 1
                    summary.statement = sprintf([summary.statement ' on %d requirement'], summary.num_requirements);
                else
                    summary.statement = sprintf([summary.statement ' on %d requirements'], summary.num_requirements);
                end
                summary.num_traces_violations = sum( any(this.traces_vals<0, 2) );
                summary.statement = sprintf([summary.statement ', %d traces have violations'], summary.num_traces_violations);
                summary.num_total_violations =  sum( sum(this.traces_vals<0) );
                if  summary.num_total_violations == 1
                    summary.statement = sprintf([summary.statement ', %d requirement violation' ], summary.num_traces_violations);
                elseif summary.num_total_violations >1
                    summary.statement = sprintf([summary.statement ', %d requirement violations total' ], summary.num_traces_violations);
                end
                
                summary.num_vacuous_sat = sum(sum(summary.requirements.rob == inf));
                
                if  summary.num_vacuous_sat == 1
                    summary.statement = sprintf([summary.statement ', %d vacuous satisfaction.' ], summary.num_vacuous_sat);
                elseif summary.num_vacuous_sat >1
                    summary.statement = sprintf([summary.statement ', %d vacuous satisfactions.' ], summary.num_vacuous_sat);
                else
                    summary.statement = [summary.statement '.'];
                end
                
                
            else
                summary.statement = [summary.statement '.'];
            end
            if isa(this.BrSet, 'BreachImportData')
                summary.file_names = this.BrSet.signalGenerators{1}.file_list;
            end
            
        end
        
        function values = GetParam(this, params, ip)
            % GetParam if not found, look into BrSet
            [idx, ifound] = FindParam(this.P, params);
            if ~exist('ip', 'var')
                ip = 1:size(this.P.pts,2);
            end
            
            if all(ifound)
                values = GetParam(this.P,idx);
                values = values(:, ip);
            else
                if ischar(params)
                    params = {params};
                end
                idx_req=  ifound~=0;
                idx_data = ifound==0;
                params_req = params(idx_req);
                params_data = params(idx_data);
                if ~isempty(params_req)
                    values(idx_req,:) = this.GetParam(params_req, ip);
                end
                for it = 1:numel(ip)
                    ipts = find(this.BrSet.P.traj_ref==ip(it),1);
                    values(idx_data,it) = this.BrSet.GetParam(params_data, ipts);
                end
            end
        end
        
        function [idx, ifound, idxB, ifoundB] = FindSignalsIdx(this, signals )
            if ischar(signals)
                signals = {signals};
            end
            idxB = zeros(1, numel(signals)); ifoundB = zeros(1, numel(signals));
            [idx, ifound] = FindSignalsIdx@BreachSet(this, signals);   %
            
            if any(~ifound)&&isa(this.BrSet, 'BreachSet')   % if not a signal of the requirement, look into BrSet (system) signals
                idx_not_found = find(~ifound);
                for isig = 1:numel(idx_not_found)
                    s0 = signals{idx_not_found(isig)};
                    aliases = this.getAliases(s0);                                                                                
                    
                    [idx_s, ifound_s] = FindParam(this.P, aliases);
                    if any(ifound_s) % one alias is the one !
                        idx(idx_not_found(isig)) = idx_s(find(ifound_s, 1));
                        ifound(idx_not_found(isig)) = 1;
                    elseif (~isempty(this.BrSet))
                        [idx_sB, ifB] = FindSignalsIdx(this.BrSet, aliases);
                        if (any(ifB))
                            idxB(idx_not_found(isig)) = idx_sB(find(ifB, 1));
                            ifoundB(idx_not_found(isig)) = 1;
                        end
                    end
                end
            end
        end
        
        function [X, idxR] = GetSignalValues(this,varargin)
            % GetSignalValues if not found, look into BrSet
            nb_traj = 0;
            X = [];
            signals = varargin{1};
            if ischar(signals)
                signals = {signals};
            end
            [idxR, ifound, idxB, ifoundB] = this.FindSignalsIdx(signals); % finds signal indexes in this or this.BrSet,
            % supports bidirectional aliases ...
            if all(ifound)
                X = GetSignalValues@BreachTraceSystem(this,idxR, varargin{2:end});
            else
                if any(~ifound&~ifoundB)
                   absent =  find(~ifound&~ifoundB,1);
                   warning('GetSignalsValues:not_found', 'Signal %s not found', signals{absent});
                end
                if any(ifound)
                    idxR = idxR(ifound==1);
                    values_req = this.GetSignalValues@BreachTraceSystem(idxR, varargin{2:end});
                    if iscell(values_req)
                        nb_traj =numel(values_req);
                    else
                        nb_traj =1;
                    end
                end
                if any(ifoundB)
                    idxB = idxB(ifoundB==1);
                    if nargin>2 && isfield(this.BrSet.P, 'traj_ref')
                        itraces = varargin{2};
                        itracesB = itraces;
                        for it= 1:numel(itraces)
                            itracesB(it) = this.P.traj{itraces(it)}.param(this.Sys.DimX+1); % uses trace_data_idx_ to recover trace in B
                        end
                        varargin{2} = itracesB;
                    end
                    
                    values_data = this.BrSet.GetSignalValues(idxB, varargin{2:end});
                    if iscell(values_data)
                        nb_traj =numel(values_data);
                    else
                        nb_traj =1;
                    end                    
                end
                
                if nb_traj>1
                    X = cell(nb_traj,1);
                    for iv = 1:nb_traj
                        if any(ifound)
                            X{iv}(ifound==1,:) = values_req{iv};
                        end
                        if any(ifoundB)
                            X{iv}(ifoundB==1,:) = values_data{iv};
                        end
                    end
                else
                    if any(ifound)
                        X(ifound==1,:) = values_req;
                    end
                    if any(ifoundB)
                        X(ifoundB==1,:) = values_data;
                    end
                end
            end
        end
         
        function dom  = GetDomain(this, param)
            if ischar(param)
                param =  {param};
            end
            [idx, found] = FindParam(this.P, param);
            for i = 1:numel(idx)
                if found(i)
                    dom(i) = this.Domains(idx(i));
                elseif ~isempty(this.BrSet)
                    [idx_BrSet, found_BrSet] = FindParam(this.BrSet.P, param(i));
                    if found_BrSet
                        dom(i) = this.BrSet.Domains(idx_BrSet);
                    else
                        dom(i) = BreachDomain();
                    end
                else
                    dom(i) = BreachDomain();
                end
            end
        end
        
        function [time, rob] = GetRobustSat(this,req,trace_num, expr)
           if nargin <= 1
               req = this.req_monitors{1};
           elseif isnumeric(req)
               req = this.req_monitors{req};
           elseif ischar(req)
               req = this.get_req_from_name(req); 
           end
           assert(isa(req, 'stl_monitor'));
                      
           if nargin <= 2
              trace_num = 1;  
           end
           traj = this.P.traj{trace_num}; 
           t = traj.time;
           idx_par_req = FindParam(this.P, req.params);
           p_in = traj.param(1, idx_par_req);
           Xin = this.GetSignalValues(req.signals_in,trace_num);           
           req.init_tXp(t,Xin,p_in);
           phi_tmp = STL_Formula('phi_tmp__', expr);
           IA_flag = false;
           [time,rob] = req.get_standard_rob(phi_tmp,t,IA_flag);
           
        end
                
        function [req, idx_req] = get_req_from_name(this, req_name)
            
            req = [];
            idx_req = 0;
            for ir = 1:numel(this.precond_monitors)
                if strcmp(this.precond_monitors{ir}.name, req_name)
                    req = this.precond_monitors{ir};
                    idx_req = -ir; 
                    return;
                end
            end
            
            for ir = 1:numel(this.req_monitors)
                if strcmp(this.req_monitors{ir}.name, req_name)
                    req = this.req_monitors{ir};
                    idx_req = ir;
                    return;
                end
            end
        end
                
        %% Display
        function varargout = disp(this)
            signals_in_st = cell2mat(cellfun(@(c) (['''' c ''', ']), this.signals_in, 'UniformOutput', false));
            signals_in_st = ['{' signals_in_st(1:end-2) '}'];
            summary = this.GetStatement();
            
            st = sprintf('BreachRequirement object for signal(s): %s. %s\n', signals_in_st, summary.statement);
            
            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end                      
            
        end
        
        function varargout = PrintFormula(this)
            st = '';
            if ~isempty(this.precond_monitors)
                st = sprintf([st '--- PRECONDITIONS ---\n']);
                for ifo = 1:numel(this.precond_monitors)
                    st = [st this.precond_monitors{ifo}.disp()];
                end
                st = [st '\n'];
            end
            st = sprintf([st '--- REQUIREMENT FORMULAS ---\n']);
            for ifo = 1:numel(this.req_monitors)
                st = [st  this.req_monitors{ifo}.disp() sprintf('\n')];
            end
            st = sprintf([st '\n']);
        
            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end
               
        end
        
        function varargout = PrintSignals(this)
            st = sprintf('---- SIGNALS IN ----\n');
            for isig = 1:numel(this.signals_in)
                sig = this.signals_in{isig};
                st = [st sprintf('%s %s\n', sig, this.get_signal_attributes_string(sig))];
            end
            st = sprintf([st '\n']);
            
            if this.P.DimX>0
                st= sprintf([st '---- SIGNALS  OUT ----\n']);
                for isig = 1:this.P.DimX
                    st = [st sprintf('%s %s\n', this.P.ParamList{isig}, this.get_signal_attributes_string(this.P.ParamList{isig}))];
                end
                st = sprintf([ st '\n']);
            end
            
            if ~isempty(this.postprocess_signal_gens)
                for iog = 1:numel(this.postprocess_signal_gens)
                    if iog==1
                        st = sprintf([st '--- POSTPROCESSING ---\n']); 
                    end
                    signals_in_st = cell2mat(cellfun(@(c) (['''' c ''', ']), this.postprocess_signal_gens{iog}.signals_in, 'UniformOutput', false));
                    signals_in_st = ['{' signals_in_st(1:end-2) '}'];
                    signals_out_st = cell2mat(cellfun(@(c) (['''' c ''', ']), this.postprocess_signal_gens{iog}.signals, 'UniformOutput', false));
                    signals_out_st = ['{' signals_out_st(1:end-2) '}'];
                    st =  [st sprintf('%s --> %s\n',signals_in_st, signals_out_st)];
                end
                st = sprintf([ st '\n']);
            end
            if ~isempty(this.sigMap)
                st = [st this.PrintAliases() '\n'];
            end
            
            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end
            
        end
        
        function varargout = PrintAll(this)
            st =this.PrintFormula();
            st = sprintf([st this.PrintSignals()]);
            st = sprintf([st this.PrintParams()]);
            
             if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
             end
             
        end
                                
        function  atts = get_signal_attributes(this, sig)
            % returns nature to be included in signature
            % should req_input, additional_test_data_signal,
            atts = {};
            
            if this.is_a_data_in(sig)
                atts =union(atts, {'data_in'});
            end
            
            [idx, found, idxB, foundB ] = this.FindSignalsIdx(sig);
            
            if foundB
                sig = this.BrSet.P.ParamList{idxB};
                atts = union(atts, this.BrSet.get_signal_attributes(sig));
                if ~this.is_a_data_in(sig)
                    atts =union(atts, {'data_other'});
                end
            elseif found
                sig = this.P.ParamList{idx};
                
                if this.is_a_requirement(sig)
                    atts =union(atts, {'quant_sat_signal'});
                end
                
                if this.is_a_violation_signal(sig)
                    atts =union(atts, {'violation_signal'});
                end
                
                if this.is_a_predicate(sig)
                    atts =union(atts, {'predicate'});
                end
                
                if this.is_a_precond_in(sig)
                    atts =union(atts, {'precond_in'});
                end
                
                if this.is_a_req_in(sig)
                    atts =union(atts, {'requirement_in'});
                end
                
                if this.is_a_postprocess_in(sig)
                    atts =union(atts, {'postprocess_in'});
                end
                
                if this.is_a_postprocess_out(sig)
                    atts =union(atts, {'postprocess_out'});
                end
            else
                atts = union(atts, {'unlinked'});
            end
        end
        
        function  atts = get_param_attributes(this, param)
            % returns nature to be included in signature
            % should req_input, additional_test_data_signal,
            atts = {};
            if ~isempty(this.BrSet)
                atts = this.BrSet.get_param_attributes(param);
            end
            
            if ismember(param, this.P.ParamList(this.P.DimX+1:end))
                atts =union(atts, {'req_param'});
            end
        end
        
        function signals = GetSignalList(this)
            if ~isempty(this.BrSet)
                signals = union(this.BrSet.P.ParamList(1:this.BrSet.P.DimX), this.P.ParamList(1:this.P.DimX), 'stable');
            else
                signals = this.P.ParamList(1:this.P.DimX);
            end
            % adds in aliases
            signals = union(signals, this.getAliases(signals), 'stable');
            
        end
        
        function SetupDiskCaching(this)
            % TODO same options as SaveResults
            
        end
        
        function params = GetParamList(this)
            if ~isempty(this.BrSet)
                params = union(this.BrSet.GetParamList(), this.P.ParamList(this.P.DimX+1:end), 'stable');
            else
                params = [this.P.ParamList(this.P.DimX+1:end)];
            end
        end
        
        function [summary, success, msg, msg_id] = SaveResults(this, folder_name, varargin)
            % Additional options
            if ~exist('folder_name', 'var')
                folder_name = '';
            else
                [success,msg,msg_id] = mkdir(folder_name);
                if ~success
                    error('Failed to create folder %s, received error id [%s], with message ''%s''', folder_name,msg_id, msg);
                end
            end
            if isempty(folder_name)
                folder_name = pwd;
            end
            options = struct('FolderName', folder_name, 'ExportToExcel', false, 'ExcelFileName', 'Results.xlsx', ...
                'IncludeSignals', [],  'IncludeParams', [], 'AddWorkflowNum',[]);
            options = varargin2struct_breach(options, varargin{:});
            
            try
                folder_name = [options.FolderName filesep this.BrSet.mdl.name '_Results_' datestr(now, 'dd_mm_yyyy_HHMM')];
            catch
                folder_name = [options.FolderName filesep 'Req_Results_' datestr(now, 'dd_mm_yyyy_HHMM')];
            end
        
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
                error('Breach:SaveResult:no_trace', 'No trace to save - run Eval command first');
                return;
            end
            
            summary = this.GetSummary(options.IncludeSignals, options.IncludeParams);
            if ~isempty(options.AddWorkflowNum);
                summary.workflow = options.AddWorkflowNum;
            end            
                        
            summary_filename = [folder_name filesep 'summary'];
            tr = this.ExportTraces(options.IncludeSignals, options.IncludeParams, 'WriteToFolder', [folder_name filesep 'traces']);
            summary.traces_list = cell(1,numel(tr));
            for it = 1:numel(tr)
                summary.traces_list{it} = tr{it}.Properties.Source;
            end
            
            i_false = find(summary.num_violations_per_trace>0);
            summary.false_traces_files = cell(1, numel(i_false));
            for it = 1:numel(i_false)
                summary.false_traces_files{it} = tr{i_false(it)}.Properties.Source;
            end
            
            save(summary_filename,'-struct', 'summary');
            
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
               if isa(this.BrSet, 'BreachSet')                   
                   if this.BrSet.AliasMap.isKey(sig)
                       aliases = union(aliases, this.BrSet.getAliases(aliases),'stable');
                   end
               end
            end
        end
        
        function varargout = PrintAliases(this)
            varargout = {};
            st='';
            if ~isempty(this.AliasMap)
                st = sprintf('---- ALIASES ----\n');
                keys = union(this.sigMap.keys(), this.sigMapInv.keys());
                printed ={};
                for ik = 1:numel(keys)
                    if this.sigMap.isKey(keys{ik})
                        sig =  this.sigMap(keys{ik});
                    else
                        sig =this.sigMapInv(keys{ik});
                    end
                    
                    [idx, found, idx_B, found_B] = this.FindSignalsIdx(sig);
                    if found
                        sig = this.P.ParamList{idx};
                        aliases = setdiff(this.getAliases(sig), sig);
                        al_st = cell2mat(cellfun(@(c) ([ c ', ']), aliases, 'UniformOutput', false));
                        al_st = al_st(1:end-2);
                    elseif found_B
                        sig = this.BrSet.P.ParamList{idx_B};
                        aliases = setdiff(this.getAliases(sig), sig);
                        al_st = cell2mat(cellfun(@(c) ([ c ', ']), aliases, 'UniformOutput', false));
                        al_st = [al_st(1:end-2) ' (model or test data)'];
                    else
                        aliases = setdiff(this.getAliases(sig), sig);
                        al_st = cell2mat(cellfun(@(c) ([ c ', ']), aliases, 'UniformOutput', false));
                        al_st = [al_st(1:end-2) ' (not linked to data)' ];
                    end
                    if ~ismember(sig, printed)
                        st =sprintf( [st sprintf('%s <--> %s\n', sig, al_st )]);
                        printed = [printed {sig} aliases];
                    end
                end
                st = sprintf([st '\n']);
            end
            
            if ~isempty(st)
                if nargout == 0
                    varargout = {};
                    fprintf(st);
                else
                    varargout{1} = st;
                end
            end
            
        end
                           
        function Concat(this,other,fast)
            if nargin<=2
                fast = false;
            end
            Pother = other.P;
            idx_trace_idx = this.P.DimX+1;
            if isfield(this.P, 'traj')&&isfield(Pother, 'traj') % param needs to be consistent with trace_idx 
                num_traj = max(this.P.traj_ref);
                for it = 1:numel(other.P.traj)
                    Pother.traj{it}.param(idx_trace_idx)=Pother.traj{it}.param(idx_trace_idx)+num_traj;
                end
                trace_idx_values = GetParam(Pother, 'data_trace_idx_');
                Pother = SetParam(Pother, 'data_trace_idx_', trace_idx_values+num_traj);                                
            end
            this.P = SConcat(this.P, Pother, fast);
            
            % wild guess:
            this.P = SetParam(this.P, this.P.DimX+1,this.P.traj_ref);
            
            this.traces_vals_precond = [this.traces_vals_precond ; other.traces_vals_precond]; 
            this.traces_vals = [this.traces_vals ; other.traces_vals]; 
            this.traces_vals_vac = [this.traces_vals_vac ; other.traces_vals_vac]; 
            
            this.val= min([ this.val other.val]);
            
        end
        
        function sigs_in = get_signals_in(this)
            sigs_in = {};
            all_sigs_in =  {};
            for ipostprocess_signal_gens = 1:numel(this.postprocess_signal_gens)
                all_sigs_in = setdiff(union(all_sigs_in,this.postprocess_signal_gens{ipostprocess_signal_gens}.signals_in, 'stable'), this.postprocess_signal_gens{ipostprocess_signal_gens}.signals, 'stable');      % remove outputs of signals generators
            end
            
            for ifo = 1:numel(this.req_monitors)
                all_sigs_in = setdiff(union(all_sigs_in,this.req_monitors{ifo}.signals_in, 'stable' ), this.req_monitors{ifo}.signals, 'stable');             % remove outputs of formula
            end
            
            for ifo = 1:numel(this.precond_monitors)
                all_sigs_in = setdiff(union(all_sigs_in, this.precond_monitors{ifo}.signals_in, 'stable'), this.precond_monitors{ifo}.signals, 'stable');             % remove outputs of formula
            end
            
            % got all sigs_in which are not sigs out
            [~, found] = this.FindSignalsIdx(all_sigs_in);
            all_sigs_in = all_sigs_in(found==0);
            
            if size(all_sigs_in,1)>1
                all_sigs_in = all_sigs_in';
            end
            
            if ~isempty(all_sigs_in)
                reps_sigs_in = all_sigs_in(1);
                aliases = this.getAliases(all_sigs_in{1});
                for is = 1:numel(all_sigs_in)
                    if ~ismember(all_sigs_in{is}, aliases)
                        reps_sigs_in= union(reps_sigs_in, all_sigs_in(is));
                        aliases = union(aliases, this.getAliases(all_sigs_in(is)));
                    end
                end
                
                for is = 1:numel(reps_sigs_in) % remove postprocess_out
                    s  = this.get_signal_attributes(reps_sigs_in{is});
                    if ~ismember('postprocess_out', s)
                        sigs_in=union(sigs_in, reps_sigs_in{is});
                    end
                end
            end
            
        end

        function names = get_monitor_names(this)
           names = cellfun(@(c)(c.name), this.req_monitors, 'UniformOutput', false); 
        end
    end
    
    
    %% Protected methods
    methods (Access=protected)
        
        function  new = getBrSet(this, varargin)
            new = true;
            switch numel(varargin)
                case 0   %  uses this
                    if ~isempty(this.BrSet)
                        B = this.BrSet;
                    else
                        error('getTraces:no_traces', 'No traces to eval requirement on.' )
                    end
                case 1  %  assumes a BreachSystem given
                    B = varargin{1};
                case 2   % time, X - easiest, sort of
                    %  orient time and X as breach usual
                    time = varargin{1};
                    X = varargin{2};
                    
                    if size(time,2)==1
                        time=time';
                    end
                    if size(time,1) ~= 1
                        error('BreachRequirement:data_inconsistent', 'time (first argument) should  be a one dimensional array.');
                    end
                    
                    if size(X,1)==numel(time)&&size(X,2)~=numel(time)
                        X = X';
                    end
                    
                    if  size(X,2) ~= numel(time)
                        error('BreachRequirement:data_inconsistent','X (second argument) should  be an array of same length as time (first arguement).');
                    end
                    
                    if numel(this.signals_in) ~= size(X,1)
                        error('BreachRequirement:data_inconsistent', 'data is inconsistent with formula' );
                    end
                    
                    traj.status = 0;
                    traj.time = varargin{1};
                    traj.X = X;
                    
                    B = BreachTraceSystem(this.signals_in,traj);
                    
                case 3 % B, params, values --> assigns values to parameters in B
                    B = varargin{1};
                    params = varargin{2};
                    values = varargin{3};
                    if ischar(params)
                        params = {params};
                    end
                    
                    % Distribute parameters to  system and requirements
                    [params_sys, i_sys] = intersect(params, B.GetParamList());
                    [params_req, i_req] = intersect(params, this.P.ParamList(this.P.DimX+1:this.P.DimP));
                    
                    if (numel(union(params_sys,params_req))) ~= numel(params) % parameter not found
                        params_not_found = setdiff(params, union(params_sys,params_req));
                        error('Parameter %s not found either as system or requirement parameter.', params_not_found{1});
                    end
                    
                    if ~isempty(params_sys)
                        B.SetParam(params_sys,values(i_sys,:),true);
                    end
                    
                    if ~isempty(params_req)
                        this.SetParam(params_req, values(i_req));
                    end
                    
            end
            
            if isa(B,'struct')   % reading one struct obtained from a SaveResult command
                B = BreachTraceSystem(B);
            end  
                                                
            % We got B, checks whether it contains necessary signals
            if isequal(B, this.BrSet)&&this.hasTraj() % same B and requirement traj computed
                new= false;
                return;
            else
                this.ResetEval(); 
            end
                            
            this.BrSet = B; 
            [ifound, ~, ifoundB, ~ ] = this.FindSignalsIdx(this.signals_in);
            if ~all(ifound|ifoundB)
               idx_not_found = find((~ifound)&(~ifoundB),1) ;
               error('BreachRequirement:Eval', 'Signal %s is not provided by data or model',this.signals_in{idx_not_found});
            end
            
            % checks parameters in B and in Req
            paramB = this.BrSet.GetParamList();
            idxR_paramR = this.P.DimX+2:this.P.DimP;
            paramR = this.P.ParamList(idxR_paramR);
            paramR_in_B = intersect(paramB, paramR);
            idxB_paramR_in_B = FindParam(B.P, paramR_in_B);
            idxR_paramR_in_B = FindParam(this.P, paramR_in_B);
            paramR_wo_B = setdiff( paramR, paramB); % parameters in R that are not in B - need to be combined
                                    
            % if this has more than one value, needs to combine them with
            % those of B
            multiple_paramR_vals  = size(this.P.pts,2)>1;
            if multiple_paramR_vals
                values_req_wo_sys =  this.GetParam(paramR_wo_B);
                B.SetParam(paramR_wo_B, values_req_wo_sys, 'combine'); % now B contains all params for requirement
                idxB_paramR = FindParam(B.P, paramR);
            end
            
            if this.eval_precond_only   % find what needs to be computed
                sigs_need = {};
                for ipr = 1:numel(this.precond_monitors)
                    sigs_need  = union(sigs_need, this.precond_monitors{ipr}.signals_in);
                end
            else
                sigs_need = this.signals_in;
            end
            
            if ~isempty(this.precond_monitors)
                for req = this.precond_monitors
                    all_inputs_req = true;
                    if ~isempty(req{1}.signals_in)
                        for sig = req{1}.signals_in
                            if ~this.is_a_model_input(sig{1})
                                all_inputs_req = false;
                                break
                            end
                        end
                    end
                    if all_inputs_req&&isa(B, 'BreachOpenSystem')
                        B.AddInputSpec(req{1}.formula_id);
                    end
                    
                end
            end
            %  checks whether we need only input signals
            all_inputs = true;
            for is =1:numel(sigs_need)
                if ~this.is_a_model_input(sigs_need{is})
                    all_inputs =false;
                    break
                end
            end
            if isa(B, 'BreachOpenSystem')
                if all_inputs
                    B.SimInputsOnly = true;
                else
                    B.SimInputsOnly = false;
                end
            end
            
            % compute traces
            if this.use_parallel
                B.SetupParallel();
            end
            
            if isa(B, 'BreachSystem')
                B.Sim();
            end
            % initialize or reset traces
            num_pts = size(B.P.pts,2);
                        
            if isfield(this.BrSet.P,'traj_ref')
                itrajs = B.P.traj_ref;
            else
                itrajs = 1:num_pts;
            end
            
            % Note: for requirements with no signal output, we could do
            % without trajR, but will see later...
            for it =1:num_pts
                trajR = struct();
                trajR.param = zeros(1,this.Sys.DimX);
                trajR.param(this.Sys.DimX+1) = itrajs(it);
                if ~multiple_paramR_vals %  parameters of R were combined previously,
                    trajR.param = [trajR.param this.P.pts(idxR_paramR,1)']; % simple case: only one pts in R otherwise...
                    trajR.param(1,idxR_paramR_in_B) = B.P.pts(idxB_paramR_in_B,it)';
                else
                    trajR.param = [trajR.param B.P.pts(idxB_paramR,it)'];
                end
                if B.hasTraj()
                    trajR.time = B.P.traj{itrajs(it)}.time;
                else
                    trajR.time = [0 1]; % some default time, not sure where else to pick a time vector
                end
                trajR.X = NaN(this.Sys.DimX, numel(trajR.time));
                this.AddTrace(trajR);
            end
        end
         
        function  [val, val_vac, traj] = pareval_req(this, req, traj, it)
        % parallel version     
            
            val_vac = NaN;            
            time = traj{it}.time;                        
            idx_sig_req = FindParam(this.P, req.signals);
            idx_par_req = FindParam(this.P, req.params);
            p_in = traj.param(1, idx_par_req)';
            if ~isempty(req.signals_in)
                Xin = this.GetSignalValues(req.signals_in, it);
            else
                Xin = [];
            end
            % checks if a signal is missing (need postprocess)
            while any(isnan(Xin))
                % should do only one iteration if postprocessing function are
                % properly ordered following dependency
                for ipp = 1:numel(this.postprocess_signal_gens)
                    psg  = this.postprocess_signal_gens{ipp};
                    Xpp_in = this.GetSignalValues(psg.signals_in, it);
                    idx_param_pp_in = FindParam(this.P, psg.params);
                    idx_Xout = FindParam(this.P, psg.signals);
                    param_pp_in = this.P.traj{it}.param(1, idx_param_pp_in);
                    [~, this.P.traj{it}.X(idx_Xout,:)]  = this.postprocessSignals(this.postprocess_signal_gens{ipp},time, Xpp_in, param_pp_in);
                end
                Xin = this.GetSignalValues(req.signals_in, it);
            end
            if ~isempty(idx_sig_req)
                [val ,traj.time, Xout, val_vac] ...
                    = this.evalRequirement(req, time, Xin, p_in);
                traj.X( idx_sig_req,:) = Xout;
            else
                [val, ~, ~, val_vac]  = this.evalRequirement(req, time, Xin, p_in);
            end
        end
        
        function  [val, val_vac] = eval_req(this, req, it)
            val_vac = NaN;            
            time = this.P.traj{it}.time;                        
            idx_sig_req = FindParam(this.P, req.signals);
            idx_par_req = FindParam(this.P, req.params);
            p_in = this.P.traj{it}.param(1, idx_par_req)';
            
            if ~isempty(req.signals_in)
                Xin = this.GetSignalValues(req.signals_in, it);
            else
                Xin = [];
            end
            % checks if a signal is missing (need postprocess)
            while any(isnan(Xin))
                % should do only one iteration if postprocessing function are
                % properly ordered following dependency
                for ipp = 1:numel(this.postprocess_signal_gens)
                    psg  = this.postprocess_signal_gens{ipp};
                    Xpp_in = this.GetSignalValues(psg.signals_in, it);
                    idx_param_pp_in = FindParam(this.P, psg.params);
                    idx_Xout = FindParam(this.P, psg.signals);
                    param_pp_in = this.P.traj{it}.param(1, idx_param_pp_in);
                    [~, this.P.traj{it}.X(idx_Xout,:)]  = this.postprocessSignals(this.postprocess_signal_gens{ipp},time, Xpp_in, param_pp_in);
                end
                Xin = this.GetSignalValues(req.signals_in, it);
            end
            if ~isempty(idx_sig_req)
                [val , this.P.traj{it}.time, Xout, val_vac] ...
                    = this.evalRequirement(req, time, Xin, p_in);
                
                this.P.traj{it}.X( idx_sig_req,:) = Xout;
            else
                [val, ~, ~, val_vac]  = this.evalRequirement(req, time, Xin, p_in);
            end
        end
        
        function  [time, Xout] = postprocessSignals(this, pp, time, Xin, pin)
            % postprocessSignals applies intermediate signals computations
            [time , Xout] = pp.computeSignals(time, Xin, pin);
        end
        
        function [val, time, Xout, val_vac] = evalRequirement(this, req, time, Xin, pin)
            % evalRequirement eval one requirement (usually one STL formula
            % monitor) on one trace, i.e., writes related signals and returns evaluations
                    
            [val, time,  Xout] = req.eval(time, Xin, pin);
            if req.set_mode('in', 'abs') 
                val_vac = req.eval(time, Xin, pin);
                req.set_mode('out', 'rel');
            else
               val_vac = NaN; 
            end
            
        end
        
        function setTraces(this, trajs_req)
            if this.hasTraj()
                this.P.traj =trajs_req;  % optimistic ... also, should fix this traj/trajs thing
            else
                for it= 1:numel(trajs_req)
                    this.AddTrace(trajs_req{it});
                end
            end
        end
        
        %% Misc
        function traj = set_signals_in_traj(this, traj, names, Xout)
            % update names with sigMap
            
            for i_sig = 1:numel(names)
                [idx1, stat] = FindParam(this.P,  names{i_sig});
                if stat == 1 % signal not found
                    traj.X(idx1, :) = Xout(i_sig,:);
                end
                
                if this.sigMap.isKey(names{i_sig})
                    idx2  = FindParam(this.P, this.sigMap(names{i_sig}));
                    traj.X(idx2, :) = Xout(i_sig,:);
                end
                
            end
            
        end
        
        function b = is_a_data_in(this, sig)
            b = false;
            aliases = this.getAliases(sig);
            b = ~isempty(intersect(aliases,this.signals_in));
        end
        
        function b = is_a_requirement(this, sig)
            b = false;
            for ifo = 1:numel(this.req_monitors)
                if isa(this.req_monitors{ifo}, 'stl_monitor')
                    if strcmp(sig, [this.req_monitors{ifo}.formula_id '_quant_sat' ])%||... 
                            %strcmp(sig, [get_id(this.req_monitors{ifo}.formula) ])
                        b = true;
                    end
                end
            end
        end
        
        function b = is_a_violation_signal(this, sig)
            b = false;
            for ifo = 1:numel(this.req_monitors)
                if isa(this.req_monitors{ifo}, 'stl_monitor')
                    
                    if strcmp(sig, [this.req_monitors{ifo}.formula_id '_violation' ])
                        b = true;
                    end
                end
            end
        end
        
        function b = is_a_predicate(this, sig)
            b = (STL_CheckID(sig) ==1);
        end
        
        function b = is_a_postprocess_out(this, sig)
            b = false;
            [idx, found] = this.FindSignalsIdx(sig);
            if found
                sig = this.P.ParamList(idx);
                for ipp = 1:numel(this.postprocess_signal_gens)
                    psg  = this.postprocess_signal_gens{ipp};
                    if ismember(sig, psg.signals)
                        b=true;
                    end
                end
            end
        end
        
        function b = is_a_postprocess_in(this, sig)
            b = false;
            [idx, found] = this.FindSignalsIdx(sig);
            if found
                sig = this.P.ParamList(idx);
                for ipp = 1:numel(this.postprocess_signal_gens)
                    psg  = this.postprocess_signal_gens{ipp};
                    if ismember(sig, psg.signals_in)
                        b=true;
                    end
                end
            end
        end
        
        function b = is_a_req_in(this, sig)
            b = false;
            [idx, found] = this.FindSignalsIdx(sig);
            if found
                sig = this.P.ParamList(idx);
                for ipp = 1:numel(this.req_monitors)
                    req  = this.req_monitors{ipp};
                    if ismember(sig, req.signals_in)
                        b=true;
                    end
                end
            end
        end
        
        function b = is_a_precond_in(this, sig)
            b = false;
            [idx, found] = this.FindSignalsIdx(sig);
            if found
                sig = this.P.ParamList(idx);
                for ipp = 1:numel(this.precond_monitors)
                    req  = this.precond_monitors{ipp};
                    if ismember(sig, req.signals_in)
                        b=true;
                    end
                end
            end
        end
        
        function b = is_a_model_input(this, sig)
            b = false;
            if ~isempty(this.BrSet)
                atts = this.BrSet.get_signal_attributes(sig);
                if ismember('model_input', atts)
                    b= true;
                end
            end
        end
        
        function b = is_a_model_output(this, sig)
            b = false;
            if ~isempty(this.BrSet)
                atts = this.BrSet.get_signal_attributes(sig);
                if ismember('model_output', atts)
                    b= true;
                end
            end
        end
        
    end
    
end