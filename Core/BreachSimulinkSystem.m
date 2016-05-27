classdef BreachSimulinkSystem < BreachOpenSystem
    % BreachSimulinkSystem Main class to interface Breach with Simulink systems
    %
    %   BrSys = BreachSimulinkSystem(mdl_name [,params, p0])
    %         
    %   Creates a BreachSystem interface to a Simulink model. 
    %
    %   Arguments: 
    %   mdl_name  -  a string naming a Simulink model.  
    %   params    -  cell array of strings | 'all'
    %   p0        -  (optional) default values for parameters
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
        lookfor_scopes = false 
    end
   
    methods
        
        function this = BreachSimulinkSystem(mdl_name, params, p0, inputfn)

            if nargin==0
                return;
            end
            
            if ~exist(mdl_name)==4  %  create Simulink system with default options
                error('BreachSimulinkSystem first argument must be the name of a Simulink model.');
            end
            
            switch nargin
                case 1,
                    this.CreateInterface(mdl_name);
                case 2,
                    this.CreateInterface(mdl_name,params);
                case 3,
                    this.CreateInterface(mdl_name, params, p0);
                case 4, 
                    this.CreateInterface(mdl_name, params, p0);
                    this.SetInputGen(inputfn);
            end
                    
            if isaSys(this.Sys) % Basically if interface was successfully created
                this.ParamRanges = [this.Sys.p(this.Sys.DimX+1:end) this.Sys.p(this.Sys.DimX+1:end)];
                this.SignalRanges = [];
                this.P = CreateParamSet(this.Sys);
                this.P.epsi(:,:) = 0;
            end
            
        end
        
        
        function CreateInterface(this, mdl, params, p0)
            
            %% Copy the model into model_breach
            
            % Give it a name
            mdl_breach = [mdl '_breach'];
            load_system(mdl);
            close_system(mdl_breach,0);
            save_system(mdl,mdl_breach);
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
            
            %% Find and Log input signals
            
            in_blks = find_system(mdl_breach,'SearchDepth',1, 'BlockType', 'Inport');
            sig_in = {};
            
            for iblk = 1:numel(in_blks)
                
                in_name = get_param(in_blks(iblk),'Name');
                in_name = regexprep(in_name,'\W','_');
                
                %Ensures port and its output line have the same name
                lh = get_param(in_blks(iblk), 'LineHandles');
                lh=lh{1}.Outport;
                set(lh, 'Name', in_name{1});
                
                %Logs input
                set(lh,'DataLoggingName', 'Use signal name', 'DataLogging',1 ,'Name', in_name{1});
                sig_in = [sig_in {in_name{1}}];
            end
            
            if isempty(sig_in)
                pu = [];
                U.params = {};
                cs.set_param('LoadExternalInput', 'off');   % Input
            else
                const_input = constant_signal_gen(sig_in);
                U.params = const_input.params;
                pu = const_input.p0';
            end
            
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
                        
            %% Scope signals
            if this.lookfor_scopes
                sig_scopes = find_scope_signals(mdl_breach);
            end
            
            %% define parameters
            exclude = {'tspan','u__','t__'};
            assignin('base','tspan', 0:1);
            if ~exist('params','var')||(strcmp('params', 'all'))          
                [params, p0] = filter_vars(mdl_breach, exclude);
            end 
            
            if ~exist('p0', 'var')||isempty(p0)
                p0 = zeros(1,numel(params));
            end
            
            params = [params U.params];
            
            %% Run the model for time 0 to check proper initialization and collect signal names
            tspan = evalin('base', 'tspan;');
            assignin('base','tspan',[0 0]);
            assignin('base','t__',0);
            assignin('base','u__',zeros(1, numel(sig_in)));
            
            simout = sim(mdl_breach);
            assignin('base','tspan',tspan);
            [~,~, signals] = simout2X(this,simout);
            
            %% Reorder sig_in
            %pos_sig_in = zeros(1, numel(sig_in));
            %for i_sig = 1:numel(sig_in)
            %    sig = sig_in{i_sig}; 
            %    pos_sig_in(i_sig) = find(strcmp(sig, signals));          
            %end
            %[~, order]= sort(pos_sig_in); 
            %sig_in = sig_in(order);
            
            %% Create the Breach structure
            p0 = [zeros(1,numel(signals)) p0 pu];
            Sys = CreateSystem(signals, params, p0'); % define signals and parameters
            
            Sys.DimU = numel(sig_in);
            Sys.InputList= sig_in;
            Sys.InputOpt = [];
            
            Sys.type= 'Simulink';
            Sys.sim = @(Sys, pts, tspan) this.sim_breach(Sys,pts, tspan);
            Sys.mdl= [mdl '_breach'];
            Sys.Dir= pwd;
            Sys.tspan = 0:t_step:t_end;
            Sys.name = Sys.mdl;  % not great..
            
            save_system(mdl_breach);
            close_system(mdl_breach);
            
            this.Sys = Sys;
            
            % Initializes InputMap and input generator
            this.InputMap = containers.Map();
            idx=0;
            for input = this.Sys.InputList
                idx = idx+1;
                this.InputMap(input{1})=idx;
            end
            
            if (~isempty(sig_in))
                InputGen = BreachSignalGen({const_input});
                this.SetInputGen(InputGen);
            end     
        end
        
        function [tout, X] = sim_breach(this, Sys, tspan, pts)
            %
            % Generic wrapper function that runs a Simulink model and collect signal
            % data in Breach format (called by ComputeTraj)
            %
            
            mdl = Sys.mdl;
            load_system(mdl);
            num_signals = Sys.DimX;
            
            params = Sys.ParamList;
            for i = 1:numel(params)-num_signals
                assignin('base',params{i+num_signals},pts(i+num_signals));
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
                simout= sim(mdl);
            catch
                s= lasterror;
                warning(['An error was returned from Simulink:' s.message '\n Returning a null trajectory']);
                tout = tspan;
                X = zeros(Sys.DimX, numel(tspan));
                return;
            end
            
            [tout, X] = simout2X(this, simout);
            
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
                            
                            % getting signal data
                            for idim =1:nbdim
                                try
                                    xdata = interp1(sig.Values.Time',double(sig.Values.Data(:,idim)),tout, 'linear','extrap');
                                    X = [X ; xdata(1,:)];
                                end
                            end
                            
                            % naming multidimensional signal= name_signal_i_
                            if nbdim==1
                                signals = {signals{:} signame};
                            else
                                for idim =1:nbdim
                                    signamei = [signame '_' num2str(idim)  '_'];
                                    signals = {signals{:} signamei};
                                end
                            end
                        end
                    end
                end
            end
        end
       
        
        function new = copy(this)
            % Instantiate new object of the same class.
            new = feval(class(this));
            
            % Copy all non-hidden properties.
            p = fieldnames(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
            if ~isempty(this.InputGenerator)
                new.InputGenerator = this.InputGenerator.copy();
                new.Sys.init_u = @(~,pts,tspan)(InitU(new,pts,tspan));
            end
            new.Sys.sim = @(Sys,pts,tspan)new.sim_breach(Sys,pts,tspan);
        end
        
        
    end
    
end
