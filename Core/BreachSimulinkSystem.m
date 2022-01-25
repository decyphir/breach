classdef BreachSimulinkSystem < BreachOpenSystem
    % BreachSimulinkSystem Main class to interface Breach with Simulink systems
    %
    %   BrSys = BreachSimulinkSystem(mdl.name [,params, p0, signals, inputfn])
    %
    %   Creates a BreachSystem interface to a Simulink model.
    %
    %   Arguments:
    %   mdl.name  -  a string naming a Simulink model.
    %   params    -  cell array of strings | 'all'
    %   p0        -  (optional) default values for parameters
    %   signals   -  specifies signals to interface
    %   inputfn   -  specifies an input generator
    %
    %   If params is not given or equal to 'all', the constructor will try
    %   to discover automatically the tunable parameters in the model.
    %   If params is empty, then the only parameters for the model will be
    %   input parameters (i.e., parameters used to generate input signals).
    %
    %   The constructor interfaces inputs, outputs and logged signals.
    %   Note that a BreachSimulinkSystem is a BreachOpenSystem, i.e., a
    %   systems that can be composed with a BreachSignalGenerator for input
    %   generation. By default, a constant input generator is created for
    %   each input of the model. Use SetInputGen method to set a different input.
    %
    
    %See also BreachOpenSystem, signal_gen
    
    properties
        FindScopes = false
        FindSignalBuilders = false
        FindFromWorkspace = false
        FindTables = false
        FindStruct = false
        MaxNumTabParam
        SimCmdArgs = {}              % argument list passed to sim command
        InputSrc                     % for each input, we match an input port or base workspace (0)
        MdlVars                      % List of variables used by the model
        SimInModelsDataFolder=false
        StopAtSimulinkError=false
        mdl
        DiskCachingRoot
    end
    
    
    methods
        
        function this = BreachSimulinkSystem(mdl_name, params, p0, signals, inputfn, varargin)
            InitBreach();
            
            if nargin ==0
                 return;                   
            end
            
            if exist('p0', 'var')&&iscell(p0)
                p0 = cell2mat(p0);
            end
            
            if ~ischar(mdl_name)
                error('BreachSimulinkSystem:wrong_argument', 'First argument of BreachSimulinkSystem must be a string naming a Simulink system.');
            end
            
            if ~exist(mdl_name)==4
                error('BreachSimulinkSystem first argument must be the name of a Simulink model.');
            end
            this.mdl.name = mdl_name;
            this.mdl.path = which(mdl_name);
            this.mdl.file_info = dir(this.mdl.path);
            this.mdl.mdl_breach_path = BreachGetModelsDataPath();
            
            %% Try reusing previous interface 
            switch nargin %try reusing previous interface
                case 1
                    BrHash = DataHash(this.mdl);
                case 2  
                    BrHash = DataHash({this.mdl, params});
                case 3  
                    BrHash = DataHash({this.mdl, params,p0});
                case 4
                    BrHash = DataHash({this.mdl, params,p0,signals});
                case 5
                    BrHash = DataHash({this.mdl, params,p0,signals,inputfn});
                otherwise
                    BrHash = DataHash([{this.mdl, params,p0, signals,inputfn} varargin]);
            end
            global BreachGlobOpt;
            if isfield(BreachGlobOpt, 'BrMap')                
                if BreachGlobOpt.BrMap.isKey(BrHash)
                    B_ = BreachGlobOpt.BrMap(BrHash);
                    this = B_.copy();
                    this.mdl.date =  datestr(now,'ddmmyy-HHMM');
                    return;
                end
            else 
                BreachGlobOpt.BrMap = containers.Map();               
            end
            
            this.mdl.date =  datestr(now,'ddmmyy-HHMM');                                                
            this.ParamSrc = containers.Map();            
            this.SetupOptions(varargin{:})
            
            switch nargin
                case 1
                    this.CreateInterface(mdl_name);
                case 2
                    this.CreateInterface(mdl_name,params);
                case 3
                    this.CreateInterface(mdl_name, params, p0);
                case 4
                    this.CreateInterface(mdl_name, params, p0, signals);
                case 5
                    this.CreateInterface(mdl_name, params, p0, signals);
                    if ~isempty(inputfn)
                        this.SetInputGen(inputfn);
                    end
                otherwise % additional options
                    this.CreateInterface(mdl_name, params, p0, signals);
                    if ~isempty(inputfn)
                        this.SetInputGen(inputfn);
                    end
            end
            
            if isaSys(this.Sys) % Basically if interface was successfully created
                this.SignalRanges = [];
                this.P = CreateParamSet(this.Sys);
                this.P.epsi(:,:) = 0;
            end
            
            %% Init domains
            this.CheckinDomain();
            
            %% Init param sources, if not done already
            for ip = this.P.DimX+1:length(this.P.ParamList)
                if  ~this.ParamSrc.isKey(this.P.ParamList{ip})
                    this.ParamSrc(this.P.ParamList{ip}) = BreachParam(this.P.ParamList{ip});
                end
            end
            
            %%
            if isempty(this.InputGenerator)
                this.InputGenerator = BreachSignalGen(constant_signal_gen({}));
            end
            
            %%            
            BreachGlobOpt.BrMap(BrHash) = this.copy();
            
        end
        
        function SetupOptions(this, varargin)
            % handles additional options for model interfacing
            options.UseDiskCaching = false;
            options.FindScopes = false;
            options.FindTables = false;
            options.FindStruct  = false;
            options.FindSignalBuilders = false;  % TODO fixme when true
            options.SimInModelsDataFolder = false;
            options.Parallel = 'off';
            options.SimCmdArgs = {};
            options.Verbose = 1;
            options.MaxNumTabParam = 10;
            options.InitFn = '';
            options = varargin2struct_breach(options, varargin{:});
            
            global BreachGlobOpt;                        
            options.DiskCachingRoot = [BreachGlobOpt.breach_dir filesep 'Ext' filesep 'ModelsData' filesep 'Cache'] ;
            
            this.UseDiskCaching = options.UseDiskCaching;
            this.DiskCachingRoot = options.DiskCachingRoot;
            this.FindScopes = options.FindScopes;
            this.FindTables = options.FindTables;
            this.FindStruct = options.FindStruct;
            this.MaxNumTabParam = options.MaxNumTabParam;
            this.SimCmdArgs = options.SimCmdArgs;
            this.verbose = options.Verbose;
            this.SimInModelsDataFolder = options.SimInModelsDataFolder;
            
            if ~isempty(options.InitFn)
                this.SetInitFn(options.InitFn);
            end
            
        end
        
        function SetupParallel(this, NumWorkers, varargin)
            % BreachSimulinkSystem.SetupParallel
            
            switch nargin
                case 1
                    this.SetupParallel@BreachSystem();
                case 2
                    this.SetupParallel@BreachSystem(NumWorkers);
                otherwise
                    this.SetupParallel@BreachSystem(NumWorkers, varargin{:});
            end
            
            pctRunOnAll('warning(''off'', ''Simulink:Commands:MdlFileChangedCloseManually'')'); % should be harmless right?
            
        end
        
        function StopParallel(this)
            if license('test','Distrib_Computing_Toolbox')
                this.StopParallel@BreachSystem();
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if ~isempty(poolobj)
                    delete(poolobj);     % not sure this is doing anything
                end
                
                this.use_parallel = 0;
                this.Sys.Parallel = 0;
            end
        end
        
        function ME=EnableFastRestart(this)
            mdl_breach = this.Sys.mdl;   
            try
                load_system(mdl_breach);
                set_param(mdl_breach, 'FastRestart', 'on');
                save_system(mdl_breach);
                ME=0;
            catch ME
               warning('BreachSimulinkSystem:no_fast_restart', ['Did not manage to enable fast restart for ' this.Sys.mdl]);
           end
        end
        
        
        %% Interface creation
        function [sig_in, sig_out, sig_fw, params, sig_build_params] = CreateInterface(this, mdl, params, p0, signals)
            %% Checks model file
            if ~isfile(mdl)
                [~,mdl] = fileparts(mdl);                
            end            
            try 
                load_system(mdl);
            catch
                error('BreachSimulinkSystem:not_a_simulink_model','%s does not exists or is not a valid Simulink filename',mdl);
            end
            
            %% Copy the model            
            %  Get Breach directory
            global BreachGlobOpt
            breach_dir = BreachGlobOpt.breach_dir;
            breach_data_dir = [breach_dir filesep 'Ext' filesep 'ModelsData' ];
            
            % Give it a name
            mdl_breach = [mdl '_breach'];            
            
            close_system(mdl_breach,0);
            save_system(mdl,[breach_data_dir filesep mdl_breach]);
            close_system(mdl,0);
            load_system(mdl_breach);
            
            %% Adjust configuration parameters of the model
            cs = getActiveConfigSet(mdl_breach);
            
            % Do not change the order of the following commands. There are dependencies between the parameters.
            cs.set_param('GenerateReport', 'off');   % Create code generation report
            cs.set_param('LaunchReport', 'off');   % Open report automatically
            cs.set_param('OptimizeBlockIOStorage', 'on');   % Signal storage reuse
            cs.set_param('ExpressionFolding', 'on');   % Eliminate superfluous local variables (expression folding)
            cs.set_param('SaveFinalState', 'off');   % Final states
            cs.set_param('SignalLogging', 'on');   % Signal logging
            cs.set_param('SaveOutput', 'on');   % Output
            
            cs.set_param('LimitDataPoints', 'off');   % Limit data points to last
            cs.set_param('LoadExternalInput', 'on');   % Input
            cs.set_param('LoadInitialState', 'off');   % Initial state
            cs.set_param('ReturnWorkspaceOutputs', 'on');   % Save simulation output as single object
            
            %%  Solver pane - times
            t_end= evalin('base',cs.get_param('StopTime'));
            if isinf(t_end)
                warning('BreachSimulinkSystem:t_end_inf', 'stop time is inf, setting to 1 instead. Use SetTime method to specify another simulation end time.');
                t_end= 1;
            end
            
            try
                min_step =str2num(cs.get_param('MinStep'));
            catch 
                min_step=[];
            end
            
            try
                max_step = str2num(cs.get_param('MaxStep'));
            catch 
                max_step=[];
            end
            
            try
                t_step= str2num(cs.get_param('FixedStep'));
            catch % default fixed step is t_end/1000, unless MaxStep is set smaller
                t_step=[];
            end

            if isempty(t_step) % makes it some default
                t_step= t_end/1000;
                if ~isempty(max_step)
                    t_step = min([t_step max_step]);
                end
                if ~isempty(min_step)
                    t_step = max([t_step min_step]);
                end
            end
            cs.set_param('StartTime', '0.0');   % Start time
            cs.set_param('StopTime', 'tspan(end)');   % Stop time
            cs.set_param('SaveTime', 'on');   % Time
            cs.set_param('TimeSaveName', 'tout');   % Time
            
            %% Data Import/Export pane
            cs.set_param('ExternalInput', '[t__, u__]');   % Input
            cs.set_param('InspectSignalLogs', 'off');   % Record and inspect simulation output
            cs.set_param('OutputSaveName', 'yout');   % Output
            cs.set_param('ReturnWorkspaceOutputsName', 'out');   % Save simulation output as single object
            cs.set_param('SaveCompleteFinalSimState', 'off');   % Save complete SimState in final state
            cs.set_param('SaveFormat', 'StructureWithTime');   % Format
            cs.set_param('SignalLoggingName', 'logsout');   % Signal logging name
            
            if (verLessThan('matlab','R2011a'))
                error('Sorry, this version of Matlab is not supported.')
            end
            cs.set_param('DSMLoggingName', 'dsmout');   % Data stores logging name
            cs.set_param('SignalLoggingSaveFormat', 'Dataset');   % Signal logging format
            
            
            %% Find and log input signals
            in_blks = find_system(mdl_breach,'SearchDepth',1, 'BlockType', 'Inport');
            nb_inputs= numel(in_blks);
            
            sig_in = cell(1, nb_inputs);
            
            for iblk = 1:nb_inputs
                
                in_name = get_param(in_blks(iblk),'Name');
                in_name = regexprep(in_name,'\W','_');
                
                %Ensures port and its output line have the same name
                lh = get_param(in_blks(iblk), 'LineHandles');
                lh=lh{1}.Outport;
                set(lh, 'Name', in_name{1});
                
                %Get port number, makes sure sig_in is in the right order
                st_port_nb = get_param(in_blks(iblk),'Port');
                port_nb = str2num(st_port_nb{1});
                
                %Logs input
                set(lh,'DataLoggingName', 'Use signal name', 'DataLogging',1 ,'Name', in_name{1});
                sig_in{port_nb} = in_name{1};
                
            end
            
            % creates default constant input generator
            if isempty(sig_in)
                pu = [];
                U.params = {};
                cs.set_param('LoadExternalInput', 'off');   % Input
                const_input = {};
            else
                const_input = constant_signal_gen(sig_in);
                U.params = const_input.params;
                pu = const_input.p0';
            end
            
            %% From workspace blocks
            sig_fw = {};
            if this.FindFromWorkspace
                ifw = find_system(mdl_breach, 'BlockType', 'FromWorkspace');
                for i = 1:numel(ifw)
                    nm = get_param(ifw{i}, 'VariableName');
                    sig_fw = {sig_fw{:}, nm};
                    % ensure consistency of signal and ToWorkspace block name -
                    % maybe not necessary, but why not?
                    line_out = get_param(ifw{i}, 'LineHandles');
                    
                    set(line_out.Outport,'Name',nm);
                    set(line_out.Outport,'DataLoggingName', 'Use signal name', 'DataLogging',1 );
                end
            end
            % creates default from_workspace input generator
            if ~isempty(sig_fw)
                fw_input = from_workspace_signal_gen(sig_fw);
            else
                fw_input = {};
            end
            default_input_gens = [const_input, fw_input];
            this.Sys.InputList= [sig_in, sig_fw]; % used by FindLoggedSignals
            this.InputSrc = zeros(1,numel(this.Sys.InputList));
            this.InputSrc(1:numel(sig_in)) = 1:numel(sig_in);
            
            %% Find outputs
            o = find_system(mdl_breach,'SearchDepth',1, 'BlockType', 'Outport');
            sig_out= {};
            for i = 1:numel(o)
                nm = regexprep(o{i},[mdl_breach '/'],'');
                nm = regexprep(nm,'\W','_');
                sig_out = {sig_out{:}, nm};
                % ensure consistency of signal and output block name
                line_out = get_param(o{i}, 'LineHandles');
                set(line_out.Inport,'Name',nm);
                set_param(o{i},'Name', nm);
            end
            
            %% To workspace blocks
            otw = find_system(mdl_breach, 'BlockType', 'ToWorkspace');
            for i = 1:numel(otw)
                nm = get_param(otw{i}, 'VariableName');
                sig_out = {sig_out{:}, nm};
                % ensure consistency of signal and ToWorkspace block name -
                % maybe not necessary, but why not?
                line_out = get_param(otw{i}, 'LineHandles');
                
                set(line_out.Inport,'Name',nm);
                set(line_out.Inport,'DataLoggingName', 'Use signal name', 'DataLogging',1 );
            end
            
            %% Scope signals
            if this.FindScopes
                sig_scopes = find_scope_signals(mdl_breach);
            end
            
            %% Signal builders
            if this.FindSignalBuilders
                % collect signal names so far
                
                % Find them
                sb_potential = find_system(mdl_breach, 'BlockType', 'SubSystem');
                sb_list={};
                for ib = 1:numel(sb_potential)
                    try
                        signalbuilder(sb_potential{ib}); % will error if not a signalbuilder - would be nice to find a better test
                        sb_list = [sb_list sb_potential{ib} ];
                    end
                end
                
                % Name them and collect parameters
                sig_build_params = {};
                sig_build_p0 = [];
                for isb = 1:numel(sb_list)
                    sb = sb_list{isb};
                    line_out = get_param(sb, 'LineHandles');
                    sb_name = get_param(sb,'Name');
                    set(line_out.Outport,'Name',sb_name);  % can a signal builder have multiple output port? do we care?
                    set(line_out.Outport,'DataLoggingName', 'Use signal name', 'DataLogging',1 );
                    sb_idx = signalbuilder(sb, 'activegroup');
                    sb_param = [sb_name '_group_idx'];
                    sig_build_params = [sig_build_params sb_param];
                    sig_build_p0(end+1) = sb_idx;
                    this.ParamSrc(sb_param) = sb;
                end
            end
            
            %% Define parameters
            exclude = {'tspan','u__','t__'};
            assignin('base','tspan', 0:1);
            
            % find all possible parameters
            if ~exist('params','var')||(exist('params','var')&&isequal(params, 'all'))
                [params, p0] = this.filter_vars(mdl_breach, exclude);
                % adds in signal_builder params
                if this.FindSignalBuilders
                    params = [params sig_build_params];
                    p0 = [p0 sig_build_p0];
                end
                
            elseif ~isempty(params)&&(~exist('p0','var') || isempty(p0))
                p0 = zeros(1,numel(params));
                
                for ip = 1:numel(params)
                    % need to check for sig_builder params
                    pname = params{ip};
                    if this.FindSignalBuilders
                        idbp = strcmp(pname, sig_build_params);
                    else
                        idbp = {};
                    end
                    
                    if ~isempty(idbp)
                        % TODO add support for signal builders in
                        % BreachParam, and a corresponding BreachParam
                        p0(ip) = 1; % default value for signal builder parameter (group idx)
                    else
                        bparam = GetSimulinkBreachParam(mdl_breach, pname);
                        this.ParamSrc(pname) = bparam;
                        p0(ip) = bparam.getValue();
                    end
                end
            elseif ~isempty(params)&&exist('p0', 'var')
                if ~isequal(size(p0), size(params))
                    error('BreachSimulinkSystem:params_p0_size_differ', 'Wrong number of default values for parameters.')
                else
                    for ip = 1:numel(params)
                        pname = params{ip};
                        bparam = GetSimulinkBreachParam(mdl_breach, pname);
                        bparam.setValue(p0(ip));
                        this.ParamSrc(pname) = bparam;
                    end
                end
            else
                params = {};
                p0=[];
            end
            
            if ~exist('p0', 'var')||isempty(p0)
                p0 = zeros(1,numel(params));
            end
            
            params = [params U.params];
            
            %% Test run model on 0 time and collect simout
            tspan = evalin('base', 'tspan;');
            if isempty(min_step)
                min_step =eps;
            end
            assignin('base','tspan',[0 51*min_step]);
            assignin('base','t__',0);
            assignin('base','u__',zeros(1, numel(this.Sys.InputList)));
            if this.SimInModelsDataFolder
                crd = pwd;
                cd(breach_data_dir);
            end
            try
                simout = sim(mdl_breach, this.SimCmdArgs{:});
            catch MException
                if this.SimInModelsDataFolder
                    cd(crd);
                end
                rethrow(MException);
            end
            if this.SimInModelsDataFolder
                cd(crd);
            end
            assignin('base','tspan',tspan);                       
            
            %% Try to turn on FastRestart
            %try
            %    set_param(mdl_breach, 'FastRestart', 'on');
            %    save_system(mdl_breach);                
            %catch ME
            %    warning('BreachSimulinkSystem:no_fast_restart', ['Did not manage to enable fast restart for ' this.Sys.mdl]);
            %end
                                    
            %% find logged signals (including inputs and outputs)
            this.Sys.mdl= mdl_breach;
            if ~exist('signals', 'var')||isempty(signals)||isequal(signals,'all')
                signals = FindLoggedSignals(this, simout);
                % Ensure inputs are at the end of signals:
                signals= setdiff(signals, sig_in);
                signals = [signals sig_in];
            else
                sig_log = FindLoggedSignals(this, simout);
                found = ismember(signals, sig_log);
                
                if ~all(found)
                    not_found = find(~found);
                    warning('BreachSimulinkSystem:signal_not_found',['Signal ' signals{not_found} ' not found in model.']);
                end
            end
            
            %% Create the Breach structure
            p0 = [zeros(1,numel(signals)) p0 pu];
            Sys = CreateSystem(signals, params, p0'); % define signals and parameters
            
            Sys.DimU = numel(sig_in);
            Sys.InputList= [sig_in sig_fw];
            Sys.InputOpt = [];
            
            Sys.type= 'Simulink';
            Sys.sim = @(Sys, pts, tspan) this.sim_breach(Sys,pts, tspan);
            Sys.mdl= [mdl '_breach'];
            Sys.Dir= breach_data_dir;
            Sys.tspan = 0:t_step:t_end;
            Sys.name = Sys.mdl;  % not great..
            
            this.Sys = Sys;
            this.P = CreateParamSet(Sys);
            
            % Initializes InputMap and input generator
            this.InputMap = containers.Map();
            idx=0;
            for input = this.Sys.InputList
                idx = idx+1;
                this.InputMap(input{1})=idx;
            end
            
            if (~isempty(default_input_gens))
                InputGen = BreachSignalGen(default_input_gens);
                this.SetInputGen(InputGen);
            end
            
            %% Setup param domains
            this.Domains = repmat(BreachDomain('double'),[1 this.Sys.DimP]);
            
            % Parameters for signalbuilder
            if this.FindSignalBuilders
                for isb = 1:numel(sb_list)
                    sb = sb_list{isb};
                    idx= FindParam(this.Sys, sig_build_params{isb});
                    [~,~,~,groupnames] = signalbuilder(sb);
                    num_groups = numel(groupnames);
                    this.Domains(idx) = BreachDomain('int',[1 num_groups]);
                end
            end
            
            %% Closing
            save_system(mdl_breach);
            
        end
        
        function sig_log = FindLoggedSignals(this, simout)
            %
            % converts a simulink output to a data structure Breach can handle
            %
            
            %% Outputs and scopes
            Vars = simout.who;
            lenVars = numel(Vars);
            sig_log = {};
            
            for iV = 1:lenVars
                Y = get(simout,Vars{iV});
                if ~isempty(Y)
                    
                    if ~strcmp(Vars{iV}, 'tout')&&~strcmp(Vars{iV},'logsout')&&(isstruct(Y))
                        for iS=1:numel(Y.signals)
                            signame = Y.signals(iS).label;
                            if ~ismember(signame,sig_log)
                                
                                nbdim = size(double(Y.signals(iS).values),2);
                                if (nbdim==1)
                                    sig_log = {sig_log{:} signame };
                                else
                                    for idim = 1:nbdim
                                        signamei = [signame '_' num2str(idim)  '_'];
                                        sig_log = {sig_log{:} signamei};
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            logs = simout.get('logsout');
            
            if ~isempty(logs)
                logs_names = logs.getElementNames();
                
                %% logs
                for ilg = 1:numel(logs_names)
                    if ~(ismember(logs_names{ilg}, sig_log))
                        signame = logs_names{ilg};
                        
                        if isempty(signame)
                            % Don't try to log signals with no name
                            continue
                        end
                        if ~ismember(signame,sig_log)
                            
                            sig = logs.getElement(signame);
                            % JOHAN CHANGE
                            try
                                if sig.numElements > 1
                                    sig = get(sig,1);

                                end
                            catch
                                % Do nothing
                            end
                            
                            try
                                nbdim = size(sig.Values.Data,2);
                            catch 
                                % Sometimes, this doesn't work
                                nbdim = length(fieldnames(sig.Values));
                            end
                            
                            % NOTE!
                            % Do we want to split multidimensional signals
                            % or not?
                            % For logged signals, currently we do NOT!
                            sig_log = {sig_log{:} signame};
                            
                            % Alternatively, if we WANT to split them, use
                            % the code below INSTEAD:
                            
                            % naming multidimensional signal= name_signal_i_
%                             if nbdim
                        end
                    end
                end
            end
        end
        
        function [vars, vals, ParamSrc] =  filter_vars(this, mdl, exclude)
            % FILTER_VAR filter variables found in Simulink models : exclude capitalized
            %   constants and change lookup tables into single variables
            %
            % Syntax: [vars vals] =  filter_var(VarsOut, exclude)
            %
            %   exclude is a regular expression patterns that should not be in the
            %   names of the variables
            %
            %   Ex: [vars vals] = filter_vars( 'model', '[A-Z]') will exclude all
            %   variable with capitalized letters in them
            %
            
            
            load_system(mdl);
            
            VarsOut = Simulink.findVars(mdl);
            mdl_ws = get_param(mdl, 'modelworkspace');
            newp = 0; % counts parameters found
            if nargin == 2
                exclude= {'tspan'};
            end
            vars ={};
            vals = [];
            for i = 1:numel(VarsOut)
                var  = VarsOut(i);
                vname = var.Name;
                if any(~cellfun(@isempty, regexp(vname, exclude)))
                    this.disp_msg(sprintf('Excluding variable %s', vname),2);
                    continue;
                end
                
                % get parameter value
                switch var.SourceType
                    case 'base workspace'
                        v = evalin('base',vname);
                        msg= sprintf('Found var %s in %s used by %s.',  var.Name, var.SourceType, var.Users{1});
                        this.disp_msg(msg,1);
                        
                    case 'model workspace'
                        v = getVariable(mdl_ws,var.Name);
                        msg= sprintf('Found var %s in %s used by %s.',  var.Name, var.SourceType, var.Users{1});
                        this.disp_msg(msg,1);
                    otherwise
                        msg= sprintf('Found var %s in a %s for block %s. This type of workspace is not yet supported.',  var.Name, var.SourceType, var.Source);
                        this.disp_msg(msg,2);
                        continue; % TODO ? support other exotic workspaces (mask workspace, etc? )
                end
                
                % is it a scalar, a struct or an array?
                if (isnumeric(v)) % scalar or array
                    new_par(v, vname);
                elseif isstruct(v)&&this.FindStruct % struct
                    fds = fieldnames(v);
                    for svi = 1:numel(fds)
                        new_par(v.(fds{svi}), [vname '.' fds{svi}]);
                    end
                end
            end
            
            function new_par(nv, nvname)
                if (isscalar(nv)) % scalar
                    newp= newp+1;
                    vars{newp}= nvname;
                    vals(newp) = nv;
                    ParamSrc{newp} = GetSimulinkBreachParam(mdl, vars{newp});
                else
                    msg = sprintf('Found %s, %d %d table - ', vname,size(v,1), size(v,2) );  % TODO verbose option
                    if this.FindTables == false
                        msg = [msg ' Set ''FindTables'' option to true to generate table parameters.'];
                        this.disp_msg(msg, 1);
                        return
                    end
                    num = 0;
                    for  i = 1:size(nv,1)
                        for j= 1:size(nv,2)
                            if (nv(i,j)~=0)&&(num<=this.MaxNumTabParam)
                                num = num+1;
                                newp=newp+1;
                                vars{newp} = [nvname '__at_' num2str(i) '_' num2str(j)];
                                vals(newp) = nv(i,j);
                                ParamSrc{newp} = GetSimulinkBreachParam(mdl, vars{newp});
                            end
                        end
                    end
                    if num>this.MaxNumTabParam
                        msg = [msg ' WARNING: number of element greater than max. Increase maximum using ''MaxNumTabParam''  option.'  ];
                    end
                    this.disp_msg(msg,1);
                end
            end
            
        end
        
        %% Simulation
        function Sim(this, tspan, U)
            
            % Refresh or setup cache if needed
            if this.UseDiskCaching
                this.SetupDiskCaching();
            end
            
            switch nargin
                case 1
                    Sim@BreachOpenSystem(this);
                case 2
                    Sim@BreachOpenSystem(this, tspan);
                case 3
                    Sim@BreachOpenSystem(this, tspan, U);
                end
            %if this.use_parallel == 0 % don't autosave in parallel mode
            %    save_system(this.Sys.mdl);
            %end
        end
        
        function [tout, X, status] = sim_breach(this, Sys, tspan, pts)
            %
            % BreachSimulinkSystem.sim_breach Generic wrapper function that runs a Simulink model and collect signal
            % data in Breach format (called by ComputeTraj)
            %
            
            status = 0; % optimistic default;
            if this.SimInModelsDataFolder
                cwd = pwd;
                if this.use_parallel
                    worker_id = get(getCurrentTask(), 'ID');
                    if ~isinteger(worker_id)
                        worker_id = 1;
                    end
                    cd([this.ParallelTempRoot filesep 'Worker' int2str(worker_id)]);
                else
                    cd(this.mdl.mdl_breach_path);
                end
            end
            % JOHAN ADDED
            % We need clear mex since otherwise, we get erroneous start
            % values for signals that need InitFunctions to be run at start
            % of simulations. 
            % We only clear mex is FastRestart is turned OFF. 
            try
                fastRestart = get_param(Sys.mdl, 'FastRestart');
            catch
                % Model does not have parameter FastRestart. Probably means
                % we are in 2013b - either way, fastRestart is considered
                % off. 
                fastRestart = 'off';
            end
            
            if strcmp(fastRestart, 'off')
                clear mex;
            end
            
            % Clear the temp Simulink storage, otherwise a huge temporary
            % file (like 100GB) will be created for many runs. 
            Simulink.sdi.clear;
            % END JOHAN ADDED
            
            mdl = Sys.mdl;
            load_system(mdl);
            num_signals = Sys.DimX;
            params = Sys.ParamList;
            
            for i = 1:numel(params)-num_signals
                pname =  params{i+num_signals};
                pval  = pts(i+num_signals);
                bparam = this.ParamSrc(pname);
                bparam.setValue(pval); % set value in the appropriate workspace
            end
            
            if ischar(tspan)
                tspan = evalin('base', tspan);
            end
            
            if isfield(Sys,'init_u')
                U = Sys.init_u(Sys.InputOpt, pts, tspan);
                assignin('base','t__',U.t);
                assignin('base', 'u__',U.u);
            end
                                    
            assignin('base','tspan',tspan);
            if numel(tspan)>2
                set_param(mdl, 'OutputTimes', 'tspan',...
                    'OutputOption','SpecifiedOutputTimes');
            else
                set_param(mdl, 'OutputTimes', 'tspan',...
                    'OutputOption','RefineOutput');
            end
            
            try
                if (this.SimInputsOnly)||(~isempty(this.InputGenerator))&&this.InputGenerator.statusMap.isKey('input_spec_false')
                    tout = this.InputGenerator.P.traj{1}.time;
                    X = NaN(Sys.DimX, numel(tout));
                    idx= this.GetInputSignalsIdx();
                    Xin = this.InputGenerator.GetSignalValues(this.Sys.InputList);
                    X(idx,:) = Xin;
                    status = -2;  % error in inputs
                else
                    simout= sim(mdl, this.SimCmdArgs{:});
                    %time_to_sim = toc;
                    %disp(['Finished simulation in ' num2str(time_to_sim) 's']);
                    [tout, X] = GetXFrom_simout(this, simout);
                end
            catch s
                if numel(tspan)>1
                    tout = tspan;
                else
                    tout = [0 tspan];
                end
                warning(['An error was returned from Simulink:' s.message '\n Returning a null trajectory']);
                disp(['WARNING: An error was returned from Simulink: ' s.message]);
                X = zeros(Sys.DimX, numel(tout));
                status =-1;
                %this.addStatus(-1, MException.identifier, MException.message);
                %if this.StopAtSimulinkError
                rethrow(s);
                %end
            end
            
            % FIXME: the following needs to be reviewed
            if ~isempty(this.InputGenerator)&&this.use_precomputed_inputs==false
                this.InputGenerator.ResetSimulations();
                this.InputGenerator.resetStatus();
            end
            
            if this.SimInModelsDataFolder
                cd(cwd);
            end
        end
        
        function [tout, X] = GetXFrom_simout(this, simout)
            %
            % converts a simulink output to a data structure Breach can handle
            %
            
            signals= this.Sys.ParamList(1:this.Sys.DimX);
            tout = simout.get('tout')';
            X=zeros(numel(signals), numel(tout));
            
            % Outputs and scopes - go over logged signals and collect those we need
            Vars = simout.who;
            lenVars = numel(Vars);
            
            for iV = 1:lenVars
                Y = get(simout,Vars{iV});
                if ~isempty(Y)
                    
                    if ~strcmp(Vars{iV}, 'tout')&&~strcmp(Vars{iV},'logsout')&&(isstruct(Y))
                        for iS=1:numel(Y.signals)
                            
                            nbdim = size(double(Y.signals(iS).values),2);
                            signame = Y.signals(iS).label;
                            if (nbdim==1)
                                [lia, loc]= ismember(signame, signals);
                                if lia
                                    xx = interp1(Y.time, double(Y.signals(iS).values),tout, 'linear','extrap') ;
                                    X(loc,:) = xx;
                                end
                            else
                                for idim = 1:nbdim
                                    signamei = [signame '_' num2str(idim)  '_'];
                                    [lia, loc]= ismember(signamei, signals);
                                    if lia
                                        xx = interp1(Y.time, double(Y.signals(iS).values(:,idim)),tout, 'linear','extrap') ;
                                        X(loc,:) = xx;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % logs - go over logged signals and collect those we need
            logs = simout.get('logsout');
            if ~isempty(logs)
                logs_names = logs.getElementNames();
                
                for ilg = 1:numel(logs_names)
                    
                    signame = logs_names{ilg};
                    sig = logs.getElement(signame);
                    % JOHAN CHANGE
                    if isempty(signame)
                        break
                    end
                     
                    try
                        if sig.numElements > 1
                            sig = get(sig,1);
                        end
                    catch
                        % Do nothing
                    end
                    
                    if ~isfield(sig.Values, 'Time') && ~isa(sig.Values, 'timeseries')
                        % sig.Values does not have a field called Time. 
                        % This means that the signal is actually a BUS
                        % which has several different signal values. 
                        % We don't handle this right now, instead we skip
                        % it and print a warning that this signal should
                        % not be used in a spec
                        disp(['BreachSimulinkSystem.m: We have tried to log the bus signal ''' sig.Name '''. Skipping logging, we cannot use bus signals in specs']);
                        continue;
                    end
                    
                    if length(sig.Values.Time) < 2
                        % We have only one element - cannot interpolate
                        % This happens e.g. for FaultModeFID_ver in
                        % CMA_CIDD, which is just a constant parameter
                    end
                    % END JOHAN CHANGE
                    nbdim = size(sig.Values.Data,2);
                    
                    if (nbdim==1)
                        [lia, loc]= ismember(signame, signals);
                        if lia
                            if length(sig.Values.Time) > 1
                                % JOHAN EDIT
                                dim = size(sig.Values.Data);
                                if length(dim) == 2
                                    % The data is 2D. We can interpolate it
                                    % the standard Breach way. 
                                    % Standard case - interpolate to fill data
                                    xx = interp1(sig.Values.Time',double(sig.Values.Data(:,1)),tout, 'linear','extrap');
                                elseif length(dim) == 3
                                    % The data is 3D. We need to figure out
                                    % which 2 dimensions to use. 
                                    [maxValue, maxDim] = max(dim);
                                    
                                    % dim is e.g. [1 1 2401]. We assert
                                    % that all dimensions OTHER than maxdim
                                    % are equal to 1. 
                                    assert(sum(dim) == maxValue + length(dim) - 1, 'All dimensions other than maxDim should be equal to 1');
                                    
                                    % We would like to look at the
                                    % dimension maxDim, as well as the
                                    % dimensions before it. To do this, we
                                    % assert that maxDim > 1. 
                                    assert(maxDim == 3, 'Need to figure out what to do if maxDim not equal to 3. Maybe we should take maxDim and the dimension AFTER it (dimension 2)? Needs specific use case');
                                    
                                    % Interpolate it in the way we know it
                                    % should work (maxDim and the
                                    % dimensions before it). 
                                    
                                    % permute() below moves the first
                                    % dimension into the thrid dimension,
                                    % essentially transforming the
                                    % dimensions to [1 2401]
                                    permutedData = permute(sig.Values.Data, [2 3 1]);
                                    xx = interp1(sig.Values.Time',double(permutedData),tout, 'linear','extrap');
                                else
                                    error('We have not defined what to do if the signal data is not 2D or 3D');
                                end
                                % END JOHAN EDIT
                            else
                                % We have only one element - cannot
                                % interpolate
                                % This happens e.g for FaultModeFID_ver in
                                % CMA_CIDD, which is just a constant
                                % parameter
                                xx = repmat(sig.Values.Data(:,1), size(tout));
                            end
                            X(loc,:) = xx;
                        end
                    else
                        for idim = 1:nbdim
                            signamei = [signame '_' num2str(idim)  '_'];
                            [lia, loc]= ismember(signamei, signals);
                            if lia
                                if length(sig.Value.Time) > 1
                                    % Standard case - interpolate to fill
                                    % data
                                    xx = interp1(sig.Values.Time', double(sig.Values.Data(:,idim)),tout, 'linear','extrap') ;
                                else
                                    % We have only one element - cannot
                                    % interpolate
                                    xx = repmat(sig.Values.Data(:,idim), size(tout));
                                end
                                X(loc,:) = xx;
                            end                  
                        end
                        
                    end
                end
            end
        end
        
        function [tout, X, signals] = simout2X(this, simout)
            %
            % converts a simulink output to a data structure Breach can handle
            %
            
            tout = simout.get('tout')';
            X=[];
            
            %% Outputs and scopes
            Vars = simout.who;
            lenVars = numel(Vars);
            signals = {};
            
            for iV = 1:lenVars
                Y = get(simout,Vars{iV});
                if ~isempty(Y)
                    
                    if ~strcmp(Vars{iV}, 'tout')&&~strcmp(Vars{iV},'logsout')&&(isstruct(Y))
                        for iS=1:numel(Y.signals)
                            signame = Y.signals(iS).label;
                            if ~ismember(signame,signals)
                                
                                nbdim = size(double(Y.signals(iS).values),2);
                                try
                                    xx = interp1(Y.time, double(Y.signals(iS).values),tout, 'linear','extrap') ;
                                catch
                                    if (nbdim==1)
                                        xx = 0*tout;
                                    else
                                        xx = zeros(numel(tout), nbdim);
                                    end
                                end
                                
                                if (nbdim==1)
                                    X = [X; xx];
                                    signals = {signals{:} signame };
                                else
                                    X = [X; xx'];
                                    for idim = 1:nbdim
                                        signamei = [signame '_' num2str(idim)  '_'];
                                        signals = {signals{:} signamei};
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            logs = simout.get('logsout');
            
            if ~isempty(logs)
                logs_names = logs.getElementNames();
                
                %% logs
                for ilg = 1:numel(logs_names)
                    if ~(ismember(logs_names{ilg}, signals))
                        signame = logs_names{ilg};
                        if ~ismember(signame,signals)
                            
                            sig = logs.getElement(signame);
                            % JOHAN CHANGE
                            try
                                if sig.numElements > 1
                                    sig = get(sig,1);
                                end
                            catch
                                % Do nothing
                            end
                            % END JOHAN CHANGE
                            nbdim = size(sig.Values.Data,2);
                            
                            % naming multidimensional signal= name_signal_i_
                            if nbdim==1
                                signals = {signals{:} signame};
                            else
                                for idim =1:nbdim
                                    signamei = [signame '_' num2str(idim)  '_'];
                                    signals = {signals{:} signamei};
                                end
                            end
                            
                            
                            % getting signal data
                            for idim =1:nbdim
                                try
                                    xdata = interp1(sig.Values.Time',double(sig.Values.Data(:,idim)),tout, 'linear','extrap');
                                    X = [X ; xdata(1,:)];
                                end
                            end
                            
                        end
                    end
                end
            end
        end
        
        function U = InitU(this,pts,tspan)
            % Computes input values
            U = InitU@BreachOpenSystem(this,pts,tspan);
            if ~isnan(U.t)
                
                idx_ports= this.InputSrc>0;
                idx_fw = find(this.InputSrc==0);
                for idx = idx_fw
                    assignin('base', this.Sys.InputList{idx}, [U.t U.u(:,idx)]);
                end
                U.u = U.u(:, idx_ports);
                
            end
        end
 
        %% Misc
        function S = GetSignature(this, varargin)
            S = GetSignature@BreachOpenSystem(this, varargin{:});
            S.mdl_info = this.mdl;
        end
        
        function OpenMdl(this)
            open_system(this.mdl.name);
        end
        
        function OpenBreachMdl(this)
            open_system(this.Sys.name);
        end
        
        function changed = get_checksum(this)
            try
                load_system(this.mdl.name);
                chs = Simulink.BlockDiagram.getChecksum(this.mdl.name);
                close_system(this.mdl.name);
                changed  = ~isequal(this.mdl.checksum, chs);
                this.mdl.checksum = chs;
            catch
                warning('BreachSimulinkSystem:get_checksum_failed', 'Simulink couldn''t compute a checksum for the model.');
                this.addStatus(0, 'BreachSimulinkSystem:get_checksum_failed', 'Simulink couldn''t compute a checksum for the model.')
            end
        end
        
        function disable_checksum_warn(this)
            warning('off', 'BreachSimulinkSystem:get_checksum_failed');
        end
        
        %% Disk Caching
        % TODO adapt to BreachSystem
        function SetupDiskCaching(this, varargin)
            %  BreachSimulinkSystem.SetupDiskCaching
            
            this.UseDiskCaching = true;
            if nargin>1
                options.DiskCachingRoot  = this.DiskCachingRoot;
                if isfield(this.Sys, 'StoreTracesOnDisk')
                    options.StoreTracesOnDisk = this.Sys.StoreTracesOnDisk;
                else
                    options.StoreTracesOnDisk = true;
                end
                options = varargin2struct_breach(options, varargin{:});
                this.DiskCachingRoot = options.DiskCachingRoot;
                this.Sys.StoreTracesOnDisk  = options.StoreTracesOnDisk;
            end
            
            % The following creates the cache folder if not done already
            this.Sys.DiskCachingFolder= this.GetCachingFolder();
            [success,~, msg_id] = mkdir(this.Sys.DiskCachingFolder);
            if success == 1
                if isequal(msg_id, 'MATLAB:MKDIR:DirectoryExists')
                    this.disp_msg(['Using existing caching folder: ' this.Sys.DiskCachingFolder],2);
                else
                    this.disp_msg(['Created caching folder:' this.Sys.DiskCachingFolder],2);
                end
            else
                this.Sys.DiskCachingFolder='';
                error(['Couldn''t create caching folder'  this.Sys.DiskCachingFolder '.']);
            end
        end
        
        function ClearDiskCache(this)
            folder = this.GetCachingFolder();
            [status, message, messageid] = rmdir(folder, 's');
            if status~=1
                if ~strcmp(messageid, 'MATLAB:RMDIR:NotADirectory')
                    error(message);
                end
            else
                this.disp_msg(['Removed cache folder '  folder],2);
            end
            if isfield(this.Sys, 'StoreTracesOnDisk')&&this.Sys.StoreTracesOnDisk
                this.ResetSimulations();
            end
        end
        
        function caching_folder_name= GetCachingFolder(this, CacheRoot)
            if nargin<=1
                CacheRoot = this.DiskCachingRoot;
            end
            mdl_hash = DataHash(this.mdl.file_info);
            caching_folder_name = [CacheRoot filesep mdl_hash];
        end
        
        function hash = get_hash(this)
            st.signals = this.GetSignalsList();
            st.params = this.GetParamList();
            
        end
            
        %% Export result
        function [summary, traces] = ExportTracesToStruct(this,i_traces, varargin)
            % BreachSimulinkSystem.ExportTracesToStruct
            
            summary = [];
            traces = [];
            if ~this.hasTraj()
                error('Breach:ExportTrace:no_trace', 'No trace to export - run Sim command first');
                return;
            end
            
            num_traces = numel(this.P.traj);
            if nargin==1||isempty(i_traces)
                i_traces = 1:num_traces;
            end
            
            % Additional options
            options = struct('FolderName', [], 'PreserveTracesOrdering', false);
            options = varargin2struct_breach(options, varargin{:});
            
            if isempty(options.FolderName)
                options.FolderName = [this.mdl.name '_Results_' datestr(now, 'dd_mm_yyyy_HHMM')];
            end
            
            %% Common stuff
            % model information
            model_info = this.mdl;
            
            % parameter names
            param_names = this.GetSysParamList();
            
            % input signal names
            signal_names= this.GetSignalNames();
            idx =  this.GetInputSignalsIdx();
            input_names = signal_names(idx);
            
            % input param names
            idxp = this.GetParamsInputIdx();
            input_params = this.P.ParamList(idxp);
            
            % signal generators
            for is  = 1:numel(input_names)
                signal_gen_types{is} = class(this.InputGenerator.GetSignalGenFromSignalName(input_names{is}));
            end
            
            % signal names
            signal_names = setdiff(signal_names, input_names);
            
            % system parameters (non-input)
            sysparams_names = setdiff(param_names, input_params);
            
            if isfield(this.P,'props_names')
                spec_names = this.P.props_names;
            end
            
            %% traces
            for it = i_traces
                
                % model info
                traces(it).model_info = model_info;
                
                % params
                traces(it).params.names = sysparams_names;
                traces(it).params.values = this.GetParam(sysparams_names,it)';
                
                % time
                traces(it).time = this.P.traj{it}.time;
                
                % input signals
                traces(it).inputs.names = input_names;
                traces(it).inputs.signal_generators = signal_gen_types;
                traces(it).inputs.params.names = input_params;
                traces(it).inputs.params.values = this.GetParam(input_params, it);
                traces(it).inputs.values =  this.GetSignalValues(input_names, it);
                
                % signals
                traces(it).outputs.names =signal_names;
                %traces(it).outputs.time = this.P.traj{it}.time;
                traces(it).outputs.values = this.GetSignalValues(signal_names, it);
                
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
            summary.model_info = model_info;
            summary.date = datestr(now);
            summary.num_traces = num_traces;
            summary.test_params.names = this.GetBoundedDomains();
            summary.input_generators = signal_gen_types;
            summary.test_params.values = this.GetParam(summary.test_params.names);
            
            summary.const_params.names = setdiff( this.P.ParamList(this.P.DimX+1:end), this.GetBoundedDomains())';
            summary.const_params.values = this.GetParam(summary.const_params.names,1)';
            
            if isfield(this.P, 'props')
                summary.specs.names = spec_names;
                if ~options.PreserveTracesOrdering
                    this.SortbyRob();
                    this.SortbySat();
                end
                summary.specs.rob = this.GetSatValues();
                summary.specs.sat = summary.specs.rob>=0;
                summary.num_sat = - sum( ~summary.specs.sat, 1  );
            end
        end
        
        function [success, msg, msg_id] = SaveResults(this, folder_name, varargin)
            % Additional options
            options = struct('FolderName', folder_name, 'SaveBreachSystem', true, 'ExportToExcel', false, 'ExcelFileName', 'Results.xlsx');
            options = varargin2struct_breach(options, varargin{:});
            
            if isempty(options.FolderName)
                folder_name = [this.mdl.name '_Results_' datestr(now, 'dd_mm_yyyy_HHMM')];
                options.FolderName = folder_name;
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
            
            [summary, traces] = this.ExportTracesToStruct([], 'PreserveTracesOrdering', options.PreserveTracesOrdering);
            %saving summary
            summary_filename = [folder_name filesep 'summary'];
            save(summary_filename,'-struct', 'summary');
            
            if  options.SaveBreachSystem
                if ischar(options.SaveBreachSystem)
                    breachsys_name = options.SaveBreachSystem;
                else
                    breachsys_name = this.whoamI;
                end
                
                breachsys_filename  = [folder_name filesep breachsys_name];
                % Need to move cache into result folder
                if this.UseDiskCaching
                    trajs = this.GetTraces();
                    this.disp_msg(['Copying cached traces to ' folder_name '\traces'], 2);
                    
                    for it = 1:numel(trajs)
                        src =  trajs{it}.Properties.Source;
                        dest = [folder_name filesep 'traces' filesep 'traj_matfile' num2str(it) '.mat'];
                        this.disp_msg([src '   --->    ' dest], 2);
                        [success,msg] = copyfile(src,dest);
                        if success==0
                            error('Copy of file %s failed with message %s', src, msg);
                        end
                        this.P.traj{it} = matfile(dest);
                    end
                end
                evalin('caller', ['save(''' breachsys_filename ''', ''' breachsys_name  ''', ''-v7.3'');'] ); % I should have written here why I'm using v7.3
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
        
        function ExportToExcel(this, excel_file)  % TODO: specialize to Simulink ? or use BreachSet.ExportToExcel
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
        
        function  st = disp(this)
            if isfield(this.P, 'traj')
                nb_traj = numel(this.P.traj);
            else
                nb_traj = 0;
            end
            [name, name_found] = this.whoamI;
            
            [~,~, ~, st_status]  = GetTraceStatus(this);
            if ~name_found
                st = ['BreachSimulinkSystem interfacing model ' this.mdl.name '. It contains ' num2str(this.GetNbParamVectors()) ' samples and ' num2str(nb_traj) ' unique traces.'];
            else
                st = ['BreachSimulinkSystem ' name ' interfacing model ' this.mdl.name '. It contains ' num2str(this.GetNbParamVectors()) ' samples and ' num2str(nb_traj) ' unique traces.'];
            end
            st = sprintf('%s %s',st, st_status);
            if nargout ==0
                fprintf('%s\n', st);
            end
        end
        
    end
    
end
