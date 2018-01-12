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
    %   p0          -  (optional) default values for parameters
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
        FindTables = false
        FindStruct = false
        MaxNumTabParam
        StoreTracesOnDisk
        SimCmdArgs = {}   % argument list passed to sim command
        InputSrc          % for each input, we match an input port or base workspace (0)
        MdlVars          % List of variables used by the model
        SimInputsOnly=false % if true, will not run Simulink model
        mdl
    end
    
    methods
        
        function this = BreachSimulinkSystem(mdl_name, params, p0, signals, inputfn, varargin)
            InitBreach();
            
            if nargin==0
                return;
            end
            
            if ~exist(mdl_name)==4  %  create Simulink system with default options
                error('BreachSimulinkSystem first argument must be the name of a Simulink model.');
            end
            this.mdl.name = mdl_name;
            this.mdl.path = which(mdl_name);
            this.mdl.date =  datestr(now,'ddmmyy-HHMM');
            this.ParamSrc = containers.Map();
            
            switch nargin
                case 1,
                    this.CreateInterface(mdl_name);
                case 2,
                    this.CreateInterface(mdl_name,params);
                case 3,
                    this.CreateInterface(mdl_name, params, p0);
                case 4,
                    this.CreateInterface(mdl_name, params, p0, signals);
                case 5
                    this.CreateInterface(mdl_name, params, p0, signals);
                    if ~isempty(inputfn)
                        this.SetInputGen(inputfn);
                    end
                otherwise, % additional options
                    this.SetupOptions(varargin{:})
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
            
        end
        
        function SetupOptions(this, varargin)
            % handles additional options for model interfacing
            options.StoreTracesOnDisk = false;
            options.FindScopes = false;
            options.FindTables = false;
            options.FindStruct  = false; 
            options.FindSignalBuilders = false;  % TODO fixme when true
            options.Parallel = 'off';
            options.SimCmdArgs = {};
            options.Verbose = 1;
            options.MaxNumTabParam = 10;
            options.InitFn = '';
            options = varargin2struct(options, varargin{:});
            
            this.StoreTracesOnDisk = options.StoreTracesOnDisk;
            this.FindScopes = options.FindScopes;
            this.FindTables = options.FindTables;
            this.FindStruct = options.FindStruct;
            this.MaxNumTabParam = options.MaxNumTabParam;
            this.SimCmdArgs = options.SimCmdArgs;
            this.verbose = options.Verbose;
            
            if ~isempty(options.InitFn)
                this.SetInitFn(options.InitFn);
            end
            
        end
        
        function SetupLogFolder(this, folder_name)
            % SetupLogFolder creates a folder to log traces
            
            mdl.checksum_hash = DataHash(this.mdl.checksum);
            if nargin<2
                st = datestr(now,'ddmmyy-HHMM');
                folder_name = [this.Sys.Dir filesep this.mdl.name '-' st];
            end
            [success,msg,msg_id] = mkdir(folder_name);
            if success == 1
                if isequal(msg_id, 'MATLAB:MKDIR:DirectoryExists')
                    this.disp_msg(['Using existing logging folder at ' folder_name]);
                else
                    this.disp_msg(['Created logging folder at ' folder_name]);
                end
                this.log_folder = folder_name;
            else
                error(['Couldn''t create folder'  folder_name '.']);
            end
            
        end
        
        function SetupParallel(this, NumWorkers, varargin)
        % BreachSimulinkSystem.SetupParallel     
        
        cluster = parcluster;
        maxNumWorkers = cluster.NumWorkers;
        
        switch nargin
            case 1
                NumWorkers = maxNumWorkers;
            case 2
                NumWorkers = min(NumWorkers, maxNumWorkers);
            otherwise
        end
        
        % check existing workers
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            currentNumWorkers = 0;
        else
            currentNumWorkers = poolobj.NumWorkers;
        end
        
        if(currentNumWorkers ~= 0 && NumWorkers ~= 0 && (currentNumWorkers ~= NumWorkers))
            this.StopParallel();
            currentNumWorkers=0; 
        end
        
        if NumWorkers~=1 && currentNumWorkers == 0
            distcomp.feature( 'LocalUseMpiexec', false );   %TODO: mathworks suggested this command. it prevents calling "Mpiexec". I'm not sure the needs of that.
            poolobj=parpool(NumWorkers); %Matlab2013b or later
        end
        this.use_parallel = NumWorkers; %Arthur comment: if Breach's original design is to let use_parallel be strictly boolean, then maybe need to create a new property to track NumWorkers
        this.Sys.Parallel = NumWorkers;
        
        % run initialization function on all workers
        if ~isempty(this.InitFn)
            pctRunOnAll(this.InitFn);
        end
        pctRunOnAll('warning(''off'', ''Simulink:Commands:MdlFileChangedCloseManually'')'); % should be harmless right?
        
        end
        
        function StopParallel(this)
            
            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if ~isempty(poolobj)
                delete(poolobj);     % not sure this is doing anything
            end

            this.use_parallel = 0;
            this.Sys.Parallel = 0;
        end
        
        %% Interface creation
        
        function [sig_in, sig_out, sig_fw, params, sig_build_params] = CreateInterface(this, mdl, params, p0, signals)
            %% Copy the model
            
            %  Get Breach directory
            global BreachGlobOpt
            breach_dir = BreachGlobOpt.breach_dir;
            breach_data_dir = [breach_dir filesep 'Ext' filesep 'ModelsData' ];
            
            % Give it a name
            mdl_breach = [mdl '_breach'];
            load_system(mdl);
            
            % Get checksum of the model
            try
                chs = Simulink.BlockDiagram.getChecksum(mdl);
                this.mdl.checksum = chs;
            catch
                warning('BreachSimulinkSystem:get_checksum_failed', 'Simulink couldn''t compute a checksum for the model.');
                this.addStatus(0, 'BreachSimulinkSystem:get_checksum_failed', 'Simulink couldn''t compute a checksum for the model.')
                this.mdl.checksum = [];
            end
            
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
            t_end= str2num(cs.get_param('StopTime'));
            try
                t_step= str2num(cs.get_param('FixedStep'));
            catch % default fixed step is t_end/1000, unless MaxStep is set smaller
                t_step=[];
            end
            
            if isempty(t_step) % makes it some default
                t_step= t_end/1000;
                try
                    maxstep = cs.get_param('MaxStep');
                    t_step = min([t_step str2num(maxstep)]);
                catch
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
                error('Sorry, this version of Matlab is too old.')
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
            ifw = find_system(mdl_breach, 'BlockType', 'FromWorkspace');
            sig_fw = {};
            for i = 1:numel(ifw)
                nm = get_param(ifw{i}, 'VariableName');
                sig_fw = {sig_fw{:}, nm};
                % ensure consistency of signal and ToWorkspace block name -
                % maybe not necessary, but why not?
                line_out = get_param(ifw{i}, 'LineHandles');
                
                set(line_out.Outport,'Name',nm);
                set(line_out.Outport,'DataLoggingName', 'Use signal name', 'DataLogging',1 );
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
            
            %% find logged signals (including inputs and outputs)
            this.Sys.mdl= mdl_breach;
            if ~exist('signals', 'var')||isempty(signals)||isequal(signals,'all')
                signals = FindLoggedSignals(this);
                % Ensure inputs are at the end of signals:
                signals= setdiff(signals, sig_in);
                signals = [signals sig_in];
            else
                sig_log = FindLoggedSignals(this);
                found = ismember(signals, sig_log);
                
                if ~all(found)
                    not_found = find(~all(found));
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
        
        function sig_log = FindLoggedSignals(this)
            %
            % converts a simulink output to a data structure Breach can handle
            %
            
            %Run the model for time 0 to check proper initialization and collect signal names
            tspan = evalin('base', 'tspan;');
            assignin('base','tspan',[0 eps]);
            assignin('base','t__',0);
            assignin('base','u__',zeros(1, numel(this.Sys.InputList)));
            
            simout = sim(this.Sys.mdl);
            assignin('base','tspan',tspan);
            
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
                        if ~ismember(signame,sig_log)
                            
                            sig = logs.getElement(signame);
                            nbdim = size(sig.Values.Data,2);
                            
                            % naming multidimensional signal= name_signal_i_
                            if nbdim==1
                                sig_log = {sig_log{:} signame};
                            else
                                for idim =1:nbdim
                                    signamei = [signame '_' num2str(idim)  '_'];
                                    sig_log = {sig_log{:} signamei};
                                end
                            end
                            
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
            switch nargin
                case 1
                    Sim@BreachOpenSystem(this);
                case 2
                    Sim@BreachOpenSystem(this, tspan);
                case 3
                    Sim@BreachOpenSystem(this, tspan, U);
            end
            if this.use_parallel == 0 % don't autosave in parallel mode
                save_system(this.Sys.mdl);
            end
        end
        
        function [tout, X] = sim_breach(this, Sys, tspan, pts)
            %
            % BreachSimulinkSystem.sim_breach Generic wrapper function that runs a Simulink model and collect signal
            % data in Breach format (called by ComputeTraj)
            %
            
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
            
            %
            % TODO: fix support for signal builder using a proper
            % BreachParam
            %keys = this.ParamSrc.keys();
            %for ik = 1:numel(keys)
            %    [ipts, p_found] = FindParam(this.Sys, keys{ik});
            %   if p_found
            %        sb = this.ParamSrc(keys{ik});
            %        signalbuilder(sb, 'activegroup', pts(ipts));
            %    end
            %end
            
            assignin('base','tspan',tspan);
            if numel(tspan)>2
                set_param(mdl, 'OutputTimes', 'tspan',...
                    'OutputOption','SpecifiedOutputTimes');
            else
                set_param(mdl, 'OutputTimes', 'tspan',...
                    'OutputOption','RefineOutput');
            end
            % save system - hopefully do nothing if nothing is to be done
            % save_system(this.Sys.mdl,[], 'OverwriteIfChangedOnDisk',true);
            
            try
                if (this.SimInputsOnly)||(~isempty(this.InputGenerator))&&this.InputGenerator.statusMap.isKey('input_spec_false')
                    tout = this.InputGenerator.P.traj{1}.time;
                    X = NaN(Sys.DimX, numel(tout));
                    idx= this.GetInputSignalsIdx();
                    Xin = this.InputGenerator.GetSignalValues(this.Sys.InputList);
                    X(idx,:) = Xin;
                else
                    simout= sim(mdl, this.SimCmdArgs{:});
                    [tout, X] = GetXFrom_simout(this, simout);
                end
            catch
                s= lasterror;
                if numel(tspan)>1
                    tout = tspan;
                else
                    tout = [0 tspan];
                end
                warning(['An error was returned from Simulink:' s.message '\n Returning a null trajectory']);
                X = zeros(Sys.DimX, numel(tout));
            end
            
            % FIXME: the following needs to be reviewed
            if ~isempty(this.InputGenerator)&&this.use_precomputed_inputs==false
                this.InputGenerator.Reset();
            end
            % save system - hopefully do nothing if nothing is to be done
            %           save_system(this.Sys.mdl,[], 'OverwriteIfChangedOnDisk',true);
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
                    nbdim = size(sig.Values.Data,2);
                    
                    if (nbdim==1)
                        [lia, loc]= ismember(signame, signals);
                        if lia
                            xx = interp1(sig.Values.Time',double(sig.Values.Data(:,1)),tout, 'linear','extrap');
                            X(loc,:) = xx;
                        end
                    else
                        for idim = 1:nbdim
                            signamei = [signame '_' num2str(idim)  '_'];
                            [lia, loc]= ismember(signamei, signals);
                            if lia
                                xx = interp1(sig.Values.Time', double(sig.Values.Data(:,idim)),tout, 'linear','extrap') ;
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
        
        
        
        %% Export result
        
        function  st = disp(this)
            if isfield(this.P, 'traj')
                nb_traj = numel(this.P.traj);
            else
                nb_traj = 0;
            end
            name = this.whoamI; 
            
            if isequal(name, '__Nobody__')
            st = ['BreachSimulinkSystem interfacing model ' this.mdl.name '. It contains ' num2str(this.GetNbParamVectors()) ' samples and ' num2str(nb_traj) ' unique traces.'];
            else
            st = ['BreachSimulinkSystem ' name ' interfacing model ' this.mdl.name '. It contains ' num2str(this.GetNbParamVectors()) ' samples and ' num2str(nb_traj) ' unique traces.'];
            end
            if nargout ==0
                disp(st);
            end
        end
        
        
        function [summary, traces] = ExportTracesToStruct(this,i_traces, varargin)
            % BreachSimulinkSystem.ExportTracesToStruct
            
            summary = [];
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
            options = varargin2struct(options, varargin{:});
            
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
                this.SortbyRob();
                this.SortbySat();
                summary.specs.rob = this.GetSatValues();
                summary.specs.sat = summary.specs.rob>=0;
                summary.num_sat = - sum( ~summary.specs.sat, 1  );
            end
        end
        
        function [success, msg, msg_id] = SaveResults(this, folder_name, varargin)
            % Additional options
            options = struct('FolderName', folder_name, 'SaveBreachSystem', true, 'ExportToExcel', false, 'ExcelFileName', 'Results.xlsx');
            options = varargin2struct(options, varargin{:});
            
            if isempty(options.FolderName)
                options.FolderName = [this.mdl.name '_Results_' datestr(now, 'dd_mm_yyyy_HHMM')];
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
    end
    
end
