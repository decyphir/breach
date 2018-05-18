classdef BreachRequirement < BreachTraceSystem
    % BreachRequirement Generic requirement class - it provides helpers to get
    % data from Breach objects so different types of constraints evaluations can easily be
    % implemented.
    
    properties
        BrSet
        ogs  % output generators
        formulas
        signals_in
        traces_vals % results for individual traces & formulas
        val            % summary evaluation for all traces & formulas
    end
    
    methods
        
        function this = BreachRequirement(formulas, ogs)
            
            %% work out arguments
            if ~iscell(formulas)
                formulas=  {formulas};
            end
            
            [signals, monitors] = get_monitors(formulas);
            if  exist('ogs', 'var')&&~isempty(ogs)
                if ~iscell(ogs)
                    ogs = {ogs};
                end
                for iogs = 1:numel(ogs)
                    signals = [signals setdiff(ogs{iogs}.signals_in, signals,'stable')];
                end
            end
            this = this@BreachTraceSystem(signals);
            
            % Add output gens
            if exist('ogs', 'var')&&~isempty(ogs)
                this.ogs = ogs;
                for iogs = 1:numel(ogs)
                    this.AddOutput(ogs{iogs});
                end
            end
            
            % Add formula monitor
            this.formulas = monitors;
            for ifo = 1:numel(monitors)
                this.AddOutput(monitors{ifo});
            end
            % Figure out what signals are required input signals
            this.ResetSigMap();
              
        end
                
        function this= SetSignalMap(this, varargin)
            % SetSignalMap maps signal names to signals needed for requirement
            % evaluation
            %
            % Input:  a map or a list of pairs or two cells
            
            this.ResetSigMap();
            
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
                            this.sigMap(varargin{1}{is}) = varargin{2}{is};
                        end
                    else
                        if ischar(varargin{1})&&ischar(varargin{2})
                            this.sigMap(varargin{1}) = varargin{2};
                        else
                            error('SetSignalMap:wrong_arg', arg_err_msg);
                        end
                    end
                otherwise
                    for is = 1:numel(varargin)/2
                        try
                            this.sigMap(varargin{2*is-1}) = varargin{2*is};
                        catch
                            error('SetSignalMap:wrong_arg', arg_err_msg);
                        end
                    end
            end
            
            this.checkSignalMap();
            this.signals_in = this.get_signals_in();
      
            if this.verbose >= 2
                this.PrintSigMap();
            end
            
        end
        
        function ResetSigMap(this)
            this.sigMap = containers.Map();
            this.signals_in = this.get_signals_in();
        end
               
        function  [val, traj] = evalTrace(this,traj)
            % evalTrace evaluation function for one trace.
            traj = this.applyOutputGens(traj);
            [val, traj] = this.getRobustSignal(traj); %  computes robustness, return at time per usual STL semantics

        end
         
        function [global_val, traces_vals] = Eval(this, varargin)
            % BreachRequirement.Eval returns evaluation of the requirement -
            % compute it for all traces available and returns min (implicit
            % conjunction)
            
            % Collect traces from context and eval them
            traces_vals = this.evalAllTraces(varargin{:});
            this.traces_vals = traces_vals;
            
            % A BreachRequirement must return a single value
            global_val = min(min(traces_vals));
            this.val = global_val;
        end
          
        function F = PlotDiagnosis(this, idx_formulas, itraj)
            if nargin<2
                idx_formulas = 1;
            end
            if nargin<3
                itraj=1;
            end
            F = BreachSignalsPlot(this, {}, itraj); % empty
            for ifo =1:numel(idx_formulas)
                this.formulas{idx_formulas(ifo)}.plot_diagnosis(F);
                title(this.formulas{idx_formulas(ifo)}.formula_id, 'Interpreter', 'None')
            end
        end
        
        function summary = GetSummary(this)
            summary.statement = sprintf('%d traces evaluated', this.CountTraces());
            summary.num_traces_evaluated =size(this.traces_vals,1);
            
            if summary.num_traces_evaluated>0
                summary.val = this.val;
                summary.traces_vals = this.traces_vals;
                summary.num_requirements = size(this.traces_vals,2);
                if summary.num_requirements == 1
                    summary.statement = sprintf([summary.statement ' on %s requirement'], summary.num_requirements);
                else
                    summary.statement = sprintf([summary.statement ' on %d requirements'], summary.num_requirements);
                end
                summary.num_traces_violations = sum( any(this.traces_vals<0, 2) );
                summary.statement = sprintf([summary.statement ', %d traces have violations'], summary.num_traces_violations);
                summary.num_total_violations =  sum( sum(this.traces_vals<0) );
                if  summary.num_total_violations == 1
                    summary.statement = sprintf([summary.statement ', %d requirement violation.' ], summary.num_traces_violations);
                elseif summary.num_total_violations >1
                    summary.statement = sprintf([summary.statement ', %d requirement violations total.' ], summary.num_traces_violations);
                else
                    summary.statement = [summary.statement '.'];
                end
            else
                summary.statement = [summary.statement '.'];
            end
            summary.num_violations_per_trace =sum(this.traces_vals<0 , 2 )';
            [~,idxm ] = sort(summary.num_violations_per_trace, 2, 'descend');
            summary.idx_traces_with_most_violations = idxm;
            
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
                values(idx_req,:) = this.GetParam(params_req, ip);
                values(idx_data,:) = this.BrSet.GetParam(params_data, ip);
            end
        end
         
        function X = GetSignalValues(this,varargin)
            % GetSignalValues if not found, look into BrSet
            signals = varargin{1};
            if ischar(signals)
                signals = {signals};
            end
            [~, ifound] = FindParam(this.P, signals);
            
            if all(ifound)
                X = this.GetSignalValues@BreachTraceSystem(varargin{:});
            else
                idx_req = ifound~=0;
                idx_data = ifound==0;
                signals_req = signals(idx_req);
                signals_data = signals(idx_data);
                if any(idx_req)
                    values_req = GetSignalValues@BreachTraceSystem(this,signals_req, varargin{2:end});
                    if iscell(values_req)
                        nb_traj =numel(values_req);
                    else
                        nb_traj =1;
                    end
                end
                if any(idx_data)
                    values_data = this.BrSet.GetSignalValues(signals_data, varargin{2:end});
                    if iscell(values_data)
                        nb_traj =numel(values_req);
                    else
                        nb_traj =1;
                    end
                end
                if nb_traj>1
                    X = cell(nb_traj,1);
                    for iv = 1:nb_traj
                        if any(idx_req)
                            X{iv}(idx_req,:) = values_req{iv};
                        end
                        if any(idx_data)
                            X{iv}(idx_data,:) = values_data{iv};
                        end
                    end
                else
                    if any(idx_req)
                        X(idx_req,:) = values_req;
                    end
                    if any(idx_data)
                        X(idx_data,:) = values_data;
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
                else
                    [idx_BrSet, found_BrSet] = FindParam(this.BrSet.P, param(i));
                    dom(i) = this.BrSet.Domains(idx_BrSet);
                end
            end
        end
        
        %% Display
        function st = disp(this)
            signals_in_st = cell2mat(cellfun(@(c) (['''' c ''', ']), this.signals_in, 'UniformOutput', false));
            signals_in_st = ['{' signals_in_st(1:end-2) '}'];
            st = sprintf('BreachRequirement object for signal(s): %s\n', signals_in_st);
            if nargout == 0
                fprintf(st);
            end
            summary = this.GetSummary();
            disp(summary.statement);
        end
        
        function PrintFormula(this)
            
           fprintf(['--- FORMULAS ---\n']);
           for ifo = 1:numel(this.formulas)
               this.formulas{ifo}.disp();
           end   
           fprintf('\n');
        end
        
        function PrintSignals(this)
            disp( '---- SIGNALS IN ----')
            for isig = 1:numel(this.signals_in)
                sig = this.signals_in{isig};
                if this.sigMap.isKey(sig)
                    fprintf('%s --> %s\n', sig , this.sigMap(sig));
                else
                    fprintf('%s\n', sig);
                end        
            end
            fprintf('\n');
           disp( '---- SIGNALS  OUT ----')
           for iog = 1:numel(this.ogs)
               signals_in_st = cell2mat(cellfun(@(c) (['''' c ''', ']), this.ogs{iog}.signals_in, 'UniformOutput', false));
               signals_in_st = ['{' signals_in_st(1:end-2) '}'];
               signals_out_st = cell2mat(cellfun(@(c) (['''' c ''', ']), this.ogs{iog}.signals, 'UniformOutput', false));
               signals_out_st = ['{' signals_out_st(1:end-2) '}'];
               fprintf('%s --> %s\n',signals_in_st, signals_out_st);
           end
   
           keys = this.sigMap.keys;
           for ik = 1:numel(keys)
               if ~ismember(keys{ik}, this.signals_in)  % internal mapping
                    fprintf('%s --> %s\n', keys{ik} , this.sigMap(keys{ik}));
               end
           end
           
        fprintf('\n');    
        end
           
        function PrintAll(this)
            this.PrintFormula();
            this.PrintSignals();
            this.PrintParams();
        end
        
        
        
        function atts = get_signal_attributes(this, sig)
            % returns nature to be included in signature
            % should req_input, additional_test_data_signal,
        atts = {}; 
        if this.is_a_requirement(sig)
            atts =[atts {'requirement'}];
        end
        if this.is_a_predicate(sig)
            atts =[atts {'predicate'}];
        end
        if this.is_a_model_input(sig)
            atts =[atts {'model_input'}];
        end
        if this.is_a_model_output(sig)
            atts =[atts {'model_output'}];
        end
        
        end
        
        function signals = GetSignalNames(this)
            if ~isempty(this.BrSet)
                signals = [this.BrSet.GetSignalNames() this.P.ParamList(1:this.P.DimX)];
            else
                signals = [this.P.ParamList(1:this.P.DimX)];
            end
        end
        
    end
    
    %% Protected methods
    methods (Access=protected)
        
        function V =evalAllTraces(this,varargin)
            % BreachRequirement.evalAllTraces collect traces and apply
            % evalTrace
            
            switch numel(varargin)
                case 0   %  uses this
                    if this.hasTraj()
                        trajs = this.P.trajs;
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
                    traj.param = this.P.pts(:,1)';
                    traj.time = varargin{1};
                    traj.X = NaN(this.Sys.DimX, numel(traj.time));
                    traj = this.set_signal_in_traj(traj, this.signals_in, X); 
                    [V, traj] = this.evalTrace(traj);
                    
                    trajs = {traj};
                
                case 3 % B, params, values --> assigns values to parameters in B
                    B = varargin{1};
                    params = varargin{2};
                    values = varargin{3}; 
                    if ischar(params)
                        params = {params};
                    end
                    
                    % Distribute parameters to  system and requirements 
                    [params_sys, i_sys] = intersect(params, B.GetParamList());
                    [params_req, i_req] = intersect(params, this.GetParamList());
                    if (numel(i_sys)+numel(i_req)) ~= numel(params) % parameter not found
                           params_not_found = setdiff(params, union(param_sys,param_req));
                           error('Parameter %s not found either as system or requirement parameter.', params_not_found{1});
                    end
                    
                    if ~isempty(params_sys)
                        B.SetParam(params_sys,values(i_sys,:));
                    end
                    B.Sim(); % FIXME - need be smarter: in particular when dealing with input constraints
                    
                    if ~isempty(params_req)
                        this.SetParam(params_req, values(i_req,:));
                    end
                    
            end
            
            
            if exist('B', 'var')
                if isa(B,'struct')   % reading one struct obtained from a SaveResult command
                    B = BreachTraceSystem(B);
                end  % here we need to handle pre-conditions, input requirements, etc
           
                this.BrSet = B;
                Xs = B.GetSignalValues(this.signals_in);
                if ~iscell(Xs)
                    Xs= {Xs};
                end
                
                % Initialize values to return
                V = zeros(numel(Xs), numel(this.formulas));
                
                % collect data necessary for formla evaluation
                for  i = 1:numel(Xs)
                    if size(this.P.pts,2) == numel(Xs)
                        trajs{i}.param = this.P.pts(:,i)';
                    elseif size(this.P.pts,2)==1
                        trajs{i}.param = this.P.pts(:,1)';
                    else
                        error('Something wrong with dimensions');
                    end
                    
                    trajs{i}.time = B.P.traj{i}.time;
                    trajs{i}.X = NaN(this.Sys.DimX, numel(trajs{i}.time));
                    trajs{i} = this.set_signal_in_traj(trajs{i}, this.signals_in, Xs{i});
                    
                    if isfield(B.P.traj{i}, 'status')
                        status = B.P.traj{i}.status;
                        
                        if status~=0
                            warning('getTraces:suspicious_status', 'Trace %d has non-zero status, indicating potentially dubious data.', i);
                        end
                    end
                    [V(i,:), trajs{i}] = this.evalTrace(trajs{i});
                end
            end
        
        
            this.setTraces(trajs);
            this.val = V;
        end
         
        function traj = applyOutputGens(this, traj)
        % applyOutputGen applies intermediate signals computations
            for iog = 1:numel(this.ogs)
                Xin = this.get_signal_from_traj(traj, this.ogs{iog}.signals_in);
                pin = traj.param(FindParam(this.P, this.ogs{iog}.params));
                [~ , Xout] = this.ogs{iog}.computeSignals(traj.time, Xin, pin);
                traj = this.set_signal_in_traj(traj, this.ogs{iog}.signals,  Xout);
            end
            
        end
        
        function [val, traj] = getRobustSignal(this,traj)
            for ifo = 1:numel(this.formulas)
                Xin = this.get_signal_from_traj(traj, this.formulas{ifo}.signals_in);
                pin = traj.param(FindParam(this.P, this.formulas{ifo}.params));
                [val(ifo),  Xout] = this.formulas{ifo}.eval(traj.time, Xin, pin);
                traj  = this.set_signal_in_traj(traj, this.formulas{ifo}.signals, Xout);
            end
        end
     
        function setTraces(this, trajs)
            if this.hasTraj()
                this.P.traj =trajs;  % optimistic ... also, should fix this traj/trajs thing
            else
                for it= 1:numel(trajs)
                    this.AddTrace(trajs{it});
                end
            end
        end
        
 %% Misc       
        function checkSignalMap(this)
            % checkSignalMap warning if a signal maps to a signal that is
            % not required by requirement
            
            values = this.sigMap.values();
            if isempty(values)
            else
                for ik = 1:numel(values)
                    idx = find(strcmp(values{ik}, this.signals_in),1);
                    if isempty(idx)
                        warning('checkSignalMap:wrong_mapping', 'Signal %s does not correspond to any signal used by requirement.', values{ik});
                    end
                end
            end
            
        end
         
        function sigs_in = get_signals_in(this)
            
            sigs_in = setdiff([this.GetSignalNames() this.sigMap.keys],  this.sigMap.values, 'stable');     % remove 'outputs' of signal maps
            
            for iogs = 1:numel(this.ogs)
                sigs_in = setdiff(sigs_in, this.ogs{iogs}.signals, 'stable');      % remove outputs of signals generators
            end
            for ifo = 1:numel(this.formulas)
            sigs_in = setdiff(sigs_in, this.formulas{ifo}.signals, 'stable');             % remove outputs of formula
            end
        end
        
        function traj = set_signal_in_traj(this, traj, names, Xout)
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
        
        function b = is_a_requirement(this, sig)
            b = false;
            for ifo = 1:numel(this.formulas)
                if strcmp(sig, get_id(this.formulas{ifo}.formula))
                    b = true;
                end
            end
        end
        
        function b = is_a_predicate(this, sig)
            b = (STL_CheckID(sig) ==1);
        end
        
        function b = is_a_model_input(this, sig)
            b = false;
            if ~isempty(this.BrSet)
                atts = this.BrSet.get_signal_attributes(sig);
                if ismember('input', atts)
                    b= true;
                end
            end
        end
        
        function b = is_a_model_output(this, sig)
            b = false;
            if ~isempty(this.BrSet)
                atts = this.BrSet.get_signal_attributes(sig);
                if ismember('output', atts)
                    b= true;
                end
            end
        end
        
    end
    
end