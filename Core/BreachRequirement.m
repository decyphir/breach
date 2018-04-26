classdef BreachRequirement < BreachTraceSystem
    % BreachRequirement Generic requirement class - it provides helpers to get
    % data from Breach objects so different types of constraints evaluations can easily be
    % implemented.
    
    properties
        ogs  % output generators
        formulas
        sigMap = containers.Map()
        signals_in
    end
    
    methods
        
        function this = BreachRequirement(formulas, ogs)
            
            if nargin==0
                return
            else
                %% work out arguments
                if ~iscell(formulas)
                    formulas=  {formulas};
                end
                signals = {};
                monitors = {};
                for ifo = 1:numel(formulas)
                    formula = formulas{ifo};
                    if isa(formula, 'char')||isa(formula, 'STL_Formula')
                        monitor = stl_monitor(formula);
                        if strcmp(get_type(monitor.formula), '=>')
                            monitor = stl_A_implies_B_monitor(formula);
                        end
                        
                    elseif isa(a, 'stl_monitor')
                        monitor = formula;
                    end
                    signals = [signals monitor.signals_in];
                    monitors = [monitors {monitor}];
                end
                
                if  exist('ogs', 'var')&&~isempty(ogs)
                    if ~iscell(ogs)
                        ogs = {ogs};
                    end
                    for iogs = 1:numel(ogs)
                        signals = [signals ogs{iogs}.signals_in ];
                    end
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
            [val, traj] = this.getRobustSignal(traj, 0); %  computes robustness, return at time per usual STL semantics 
            
        end
         
        function [global_val, trace_vals] = Eval(this, varargin)
            % BreachRequirement.Eval returns evaluation of the requirement -
            % compute it for all traces available and returns min (implicit
            % conjunction)
            
            % Collect traces from context and eval them
            trace_vals = this.evalAllTraces(varargin{:});
            
            % A BreachRequirement must return a single value
            global_val = min(min(trace_vals));
            
        end
          
        function PlotDiagnosis(this)
        % Proof of concept version
            traj = this.P.traj{1};
            Xin = this.get_signal_from_traj(traj, this.formulas{1}.signals_in);
            pin = traj.param(FindParam(this.P, this.formulas{1}.params));
        
            num_phi = numel(this.formulas);
            for ifo = 1:num_phi
                subplot(num_phi,1,ifo)
                hold on;
                this.PlotSignals(this.signals_in,[], 'b',true);  % on same axis
                ax(ifo) = this.formulas{ifo}.plot_diagnosis(traj.time, Xin, pin, 'compact');
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
        end
        
        function PrintFormula(this)
            
            st = sprintf(['--- FORMULAS ---\n']);
           for ifo = 1:numel(this.formulas)
            phi = this.formulas{ifo}.formula;
            st = sprintf([st  '%s := %s'], get_id(phi), disp(phi,1));
            
            if ~strcmp(get_type(phi),'predicate')
                st = [st '\n  where\n'];
                predicates = STL_ExtractPredicates(phi);
                for ip = 1:numel(predicates)
                    st =   sprintf([ st '%s := %s \n' ], get_id(predicates(ip)), disp(predicates(ip)));
                end
                st = [st '\n'];
            end
           end
           
            if nargout ==0
                fprintf(st);
            end
            
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
        
        function PrintSigMap(this)
            st = 'Signals Map:\n';
            
            keys  = this.sigMap.keys;
            for ip = 1:numel(keys)
                st =   sprintf([ st '%s ---> %s\n' ],keys{ip}, this.sigMap(keys{ip}));
            end
            fprintf(st);
            
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
                    
            end
            
            if exist('B', 'var')
                if isa(B,'struct')
                    if isfield(B,'names')
                        % get time
                        idx_time = find(strcmpi(B.names, 'TIME'),1);
                        time = B.values(idx_time,:);
                        
                        % get signal values
                        Xs = zeros(numel(this.signals_in),numel(time));
                        for isig = 1:numel(this.signals_in)
                            idx_sig = find(strcmp(B.names, this.signals_in{isig}),1);
                            Xs(isig,:) = B.values(idx_sig,:);
                        end
                        
                        traj.status = 0;
                        traj.param = this.P.pts(:,1)';
                        traj.time = time;
                        traj.X = NaN(this.Sys.DimX, numel(traj.time));
                        traj = this.set_signal_in_traj(traj, this.signals_in, Xs);
                        [V, traj] = this.evalTrace(traj);
                        trajs = {traj};
                        
                    end
                else
                    Xs = B.GetSignalValues(this.signals_in);
                    if ~iscell(Xs)
                        Xs= {Xs};
                    end
                    
                    % Initialize values to return
                    V = zeros(1,numel(Xs));
                    
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
                        [V(i), trajs{i}] = this.evalTrace(trajs{i});
                    end
                end
            end
            
            this.setTraces(trajs);
            
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
        
        function [val, traj] = getRobustSignal(this,traj, tau)
            for ifo = 1:numel(this.formulas)
                Xin = this.get_signal_from_traj(traj, this.formulas{ifo}.signals_in);
                pin = traj.param(FindParam(this.P, this.formulas{ifo}.params));
                [~,  Xout] = this.formulas{ifo}.computeSignals(traj.time, Xin, pin);
                traj  = this.set_signal_in_traj(traj, this.formulas{ifo}.signals, Xout);
                val(ifo) = interp1(traj.time,Xout, tau, 'previous');
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
    end
    
end