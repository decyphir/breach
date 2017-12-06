classdef BreachOpenSystem < BreachSystem
    % BreachOpenSystem  a BreachSystem derivated class with an input generator.
    %
    %   BreachOpenSystem Properties
    %        InputGenerator - BreachSystem generating inputs.
    %
    %   BreachOpenSystem Methods
    %        SetInputGen - takes a BreachSystem as argument and makes it the input generator.
    %                      Can be seen as serial composition of two BreachSystems.
    %        Sim         - the Sim method for BreachOpenSystems accepts
    %                      input as a third argument, when given, it bypasses
    %                      the input generator. The input format is an array
    %                      where the first column is time.
    %
    %See also signal_gen
    
    properties
        InputMap       % Maps input signals to idx in the input generator
        InputGenerator % BreachSystem responsible for generating inputs
        use_precomputed_inputs = false % if true, will fetch traces in InputGenerator  
    end
    
    methods
        
        function  this = BreachOpenSystem(varargin)
            this=this@BreachSystem(varargin{:});
        end
        
        function Sim(this,tspan,U)

            evalin('base', this.InitFn);
            this.CheckinDomainParam();
            
            if this.use_precomputed_inputs % Input generator drives the thing 
                ig_params  = this.InputGenerator.GetParamsSysList();
                all_pts_u = this.InputGenerator.GetParam(ig_params);
                this.SetParam(ig_params, all_pts_u);
            end
            
            if ~exist('tspan','var')
                tspan = this.Sys.tspan;
            end
            
            % Checks wether we're logging to folder
            do_compute = 1; 
            if ~isempty(this.log_folder)
               % See if this simulation was performed already
               hash = DataHash({this.P.ParamList, this.P.pts, tspan});
               log_filename = [this.log_folder filesep 'Br' hash '.mat'];
               
               if exist(log_filename, 'file')
                  % additional check? 
                  load(log_filename);
                  this.P = Br.P;
                  this.disp_msg('Reading trace from log file.', 2);
                  do_compute=0;
               else
                  % TODO figure out for each trace if a log file exist 
               
               end                 
            end
            
            if (do_compute)
            Sys = this.Sys;
            if exist('U','var') % in this case, the InputGenerator becomes a trace object
                % TODO: handles multiple input signals - or use an
                % from_workspace_signal_gen?
                
                if isnumeric(U)
                    DimU = this.InputMap.Count();
                    if size(U, 2)~=DimU+1;
                        err_msg= fprintf('Input must be an array with %d columns, first one being time.',DimU);
                        error(err_msg);
                    end
                    Us.t = U(:,1);
                    Us.u = U(:,2:end);
                else
                    Us = U;
                end
                InputGen = BreachTraceSystem(this.InputMap.keys,U);
                this.SetInputGen(InputGen);
                Sys = this.Sys;
                Sys.init_u = @(~, pts, tspan) (Us);
            end
            
            this.P = ComputeTraj(Sys, this.P, tspan);
            this.CheckinDomainTraj();
            
            % Write log file 
            if ~isempty(this.log_folder)
               this.disp_msg('Writing to log file.', 2);
               
               for itraj = 1:numel(this.P.traj)
                   hash_traj = DataHash({this.Sys.ParamList, this.P.traj{itraj}.param, this.P.traj{itraj}.time});
                   log_traj_filename = [this.log_folder filesep 'traj_' hash_traj '.mat'];
                   log_traj = matfile(log_traj_filename);
                   log_traj.param = this.P.traj{itraj}.param;
                   log_traj.time = this.P.traj{itraj}.time;
                   log_traj.X = this.P.traj{itraj}.X;
                   this.P.traj{itraj} = log_traj;
                   this.P.traj{itraj}.Properties.Writable= false;
               end
               Br = this.copy();
               save(log_filename,'Br');
            end
            
            this.CheckinDomainTraj();
            end
        end
        
        % we merge parameters of the input generator with those of the
        % system, but keep both BreachObjects
        function SetInputGen(this, IG)
            % BreachOpenSystem.SetInputGen Attach a BreachSystem as input generator.
                        
            % look for property parameters and save them
            PropParams={};
            if ~isempty(this.P)
                PropParams = this.P.ParamList(this.P.DimP+1:end);
                PropParamsValues = GetParam(this.P, PropParams);
            end
            
            inputs = this.Sys.InputList;
            
            % First remove parameters from previous input generator
            Sys = this.Sys;
            idx_prev_inputs = this.GetParamsInputIdx();
            if (~isempty(idx_prev_inputs))
                idx_not_prev_inputs = boolean(ones(1, numel(Sys.ParamList)));
                idx_not_prev_inputs(idx_prev_inputs) = 0;
                Sys.ParamList = Sys.ParamList(idx_not_prev_inputs);
                Sys.p = Sys.p(idx_not_prev_inputs);
                Sys.DimP = numel(Sys.ParamList);
                this.Sys = Sys;
                this.Domains = this.Domains(idx_not_prev_inputs); % remove domains as well
            end
            
            if ischar(IG)
                pref = 'UniStep';
                if regexp(IG, [pref '[0-9]+'])
                    cp = str2num(IG(numel(pref)+1:end));
                    IG = struct('type','UniStep','cp', cp*ones(1, numel(inputs)));
                else
                    pref = 'VarStep';
                    if regexp(IG, [pref '[0-9]+'])
                        cp = str2num(IG(numel(pref)+1:end));
                        IG = struct('type','VarStep','cp', cp*ones(1, numel(inputs)));
                    end
                end
            end
            
            % IG can be a struct, a signal generator, or a BreachSystem
            DimU = this.InputMap.Count;
            if (isstruct(IG))
                if ~isfield(IG,'type')
                    error('Input generator must be a struct with fields ''type'', ''cp'' at least')
                end
                
                if isscalar(IG.cp)
                    IG.cp = IG.cp*ones(1, DimU);
                end
                
                if ~isfield(IG,'method')
                    IG.method= 'previous';
                end
                
                if ischar(IG.method)
                    IG.method = {IG.method};
                end
                
                if numel(IG.method)==1
                    method = IG.method{1};
                    IG.method = cell(1,DimU);
                    for iu = 1:DimU
                        IG.method{iu} = method;
                    end
                end
                
                switch(IG.type)
                    case 'UniStep'
                        sg = fixed_cp_signal_gen(inputs, IG.cp, IG.method);
                        IG = BreachSignalGen({sg});
                    case 'VarStep'
                        sg = var_cp_signal_gen(inputs, IG.cp, IG.method);
                        IG = BreachSignalGen({sg});
                end
            elseif iscell(IG)
                mm = methods(IG{1});
                if any(strcmp(mm,'computeSignals'))
                    IG = BreachSignalGen(IG);
                else
                    if ~any(strcmp(mm, 'Sim'))
                        error('Input generator should be a struct, a signal_gen, a cell array of signal_gen, or a BreachSignalGen object');
                    end
                end
            else
                mm = methods(IG);
                if any(strcmp(mm,'computeSignals'))
                    IG = BreachSignalGen({IG});
                else
                    if ~any(strcmp(mm, 'Sim'))
                        error('Input generator should be a struct, a signal_gen derived object or a BreachSignalGen object');
                    end
                end
            end
            
            % Check Consistency - IG must construct signals for all signals in this.InputList
            for input = this.InputMap.keys
                idx= FindParam(IG.Sys,input);
                if idx<=IG.Sys.DimX
                    this.InputMap(input{1}) = idx;
                else
                    error(['Input ' input{1} ' is not provided by input generator.']);
                end
            end
            
            this.InputGenerator = IG;
            
            % Adds parameters for new input generator
            i_params = IG.Sys.DimX+1:IG.Sys.DimP;
            this.Sys = SetParam(this.Sys, IG.Sys.ParamList(i_params), IG.P.pts(i_params,1));
            this.Sys.DimP = numel(this.Sys.ParamList);
            
            % Resets P and ranges
            this.SignalRanges = [];
            this.P = CreateParamSet(this.Sys);
            this.P.epsi(:,:) = 0;
            
            % Restore property parameter
            if ~isempty(PropParams)
                this.P = SetParam(this.P, PropParams, PropParamsValues);
            end
            % Sets the new input function for ComputeTraj
            % FIXME?: tilde?
            this.Sys.init_u = @(~, pts, tspan) (InitU(this,pts,tspan));
            
            % Copy Domains of input generator
            IGdomains = IG.GetDomain(IG.P.ParamList);
            for ip = 1:numel(IG.P.ParamList)
                parami = IG.P.ParamList{ip};
                [~, found] = FindParam(this.P, parami);
                if (found)
                    this.SetDomain(parami, IGdomains(ip));
                end
            end
            
            % Copy or init ParamSrc
            
            %% Init param sources, if not done already
            for ip = IG.P.DimX+1:numel(IG.P.ParamList)
                p  =IG.P.ParamList{ip};
                if IG.ParamSrc.isKey(p)
                    this.ParamSrc(p) = IG.ParamSrc(p);
                end
                if  ~this.ParamSrc.isKey(p)
                   this.ParamSrc(p) = BreachParam(p); 
                end
            end

        end
        
        function [params, idx] = GetParamsSysList(this)
            idx_inputs = this.GetParamsInputIdx();
            if isempty(idx_inputs)
                idx = this.Sys.DimX+1:this.Sys.DimP;
            else
                idx = this.Sys.DimX+1:idx_inputs(1)-1;
            end
            params = this.Sys.ParamList(idx);
        end
            
        function idx = GetParamsInputIdx(this)
            if isempty(this.InputGenerator)
                [~, idx] = FindParamsInput(this.Sys);
            else
                ig_params = this.InputGenerator.Sys.DimX+1:this.InputGenerator.Sys.DimP;
                idx = FindParam(this.Sys, this.InputGenerator.Sys.ParamList(ig_params));
            end
        end
        
        function AddInputSpec(this, varargin)
            this.InputGenerator.AddSpec(varargin{:});
        end
        
        function SetInputSpec(this, varargin)
            this.InputGenerator.SetSpec(varargin{:});
            this.InputGenerator.AddSpec(varargin{:});
        end
        
        % calling the Input generator -
        % there might be saving to do if inputs are pre-generated
        function U = InitU(this, pts, tspan)
            
            idx_u = this.GetParamsInputIdx();
            pts_u = pts(idx_u);
            [~, ig_params]  = this.InputGenerator.GetParamsSysList();
            
            if this.use_precomputed_inputs
                all_pts_u = this.InputGenerator.GetParam(ig_params);
                for ip =1:size(all_pts_u, 2)
                    if isequal(pts_u, all_pts_u(:,ip))
                        i_traj  = ip;
                        break
                    end
                end
                U.t = this.InputGenerator.P.traj{i_traj}.time';
                U.u= this.InputGenerator.P.traj{i_traj}.X';  % assumes only input generator has only the input signals
                                                                              % I know this will break on me eventually, hopefully not today   
            else
                
                this.InputGenerator.P = SetParam(this.InputGenerator.P,ig_params,pts_u); % FIXME: legacy use
                this.InputGenerator.Sim(tspan);
                
                % if an inputspec is violated, sabotage U into NaN
                if ~isempty(this.InputGenerator.Specs)
                    rob = this.InputGenerator.CheckSpec();
                    if rob<0
                        this.InputGenerator.addStatus(1,'input_spec_false', 'A specification on inputs is not satisfied.')
                        U.t=NaN;
                        U.u=NaN;
                        return;
                    end
                end
                
                U.t = this.InputGenerator.P.traj{1}.time;
                if size(U.t, 1)==1
                    U.t=U.t';
                end
                
                U.u = zeros(numel(U.t), this.InputGenerator.Sys.DimX);
                idx_mdl = 0;
                for input= this.Sys.InputList % Sys.InputList is in the same order as the model
                    idx_mdl = idx_mdl+1;
                    idx =  FindParam(this.InputGenerator.P, input{1});
                    U.u(:,idx_mdl) = this.InputGenerator.P.traj{1}.X(idx,:)';
                end
            end
        end
        
        % not sure why I need a special Concat operator here.
        function this = Concat(this,other)
            
            if isa(this.InputGenerator, 'BreachTraceSystem')
                % TODO Concat other with more than one trace...
                trace = [other.InputGenerator.P.traj{1}.time' other.InputGenerator.P.traj{1}.X'];
                this.InputGenerator.AddTrace(trace);
                % Using SetParam here erases the trajectory...
                i_trace_id = FindParam(other.P, 'trace_id');
                other.P.pts(i_trace_id,1) = numel(this.P.traj)+1;
                other.P.traj{1}.param(i_trace_id) = numel(this.P.traj)+1;
                this.P = SConcat(this.P, other.P);
            elseif isa(this.InputGenerator, 'BreachSystem')
                this.InputGenerator.P= SConcat(this.InputGenerator.P, other.InputGenerator.P);
                this.P = SConcat(this.P,other.P);
            else % ignore InputGenerator .. 
                this.P = SConcat(this.P,other.P);
            end
            
        end
        
        function hsi = SetInputGenGUI(this)
            hsi=  signal_gen_gui(this);
        end
        
        function idx = GetInputSignalsIdx(this)
            idx = FindParam(this.Sys, this.Sys.InputList);
        end
        
        function PrintSignals(this)
            if isempty(this.SignalRanges)
                disp( 'Signals:')
                disp( '-------')
                for isig = 1:this.Sys.DimX
                    if any(strcmp(this.Sys.ParamList{isig}, this.Sys.InputList))
                        fprintf('%s %s (Input)\n', this.Sys.ParamList{isig}, this.Domains(isig).short_disp(1));
                    else
                        fprintf('%s %s\n', this.Sys.ParamList{isig}, this.Domains(isig).short_disp(1));
                    end
                end
            else
                
                fprintf('Signals (in range estimated over %d simulations):\n', numel(this.P.traj))
                disp('-------')
                for isig = 1:this.Sys.DimX-this.Sys.DimU
                    fprintf('%s in  [%g, %g]\n', this.Sys.ParamList{isig}, this.SignalRanges(isig,1),this.SignalRanges(isig,2));
                end
                for isig =  this.Sys.DimX-this.Sys.DimU+1:this.Sys.DimX
                    fprintf('%s (Input) in  [%g, %g]\n', this.Sys.ParamList{isig}, this.SignalRanges(isig,1),this.SignalRanges(isig,2));
                end
            end
            disp(' ')
        end
        
    end
    
end


