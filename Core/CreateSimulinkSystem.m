function Sys = CreateSimulinkSystem(mdl, signals, params, p0,  inputfn, pu0, options )
% CreateSimulinkSystem  Create a Breach system structure from a simulink model
%
% Synopsis: Sys = CreateSimulinkSystem(mdl, signals, params, p0, inputfn, pu0)
%
%  - mdl       name of the Simulink model - note that breach will use its
%              own copy of this model
%  - signals   can be either 'logged' (default) or {'s1','s2', ..., 'sn'}
%              where s1, s2,..., sn are logged signals in the model
%              if signals=='logged', Breach will find and use all logged
%              signals in the model. In either case, Breach will look for
%              and monitor input and output signals
%  - params    {'p1',...,'pn'} where pi are tunable parameters in mdl. If
%              empty or absent, Breach will try to find tunable parameters
%  - p0        default values for the tunable parameters (0 if empty)
%  - inputfn   input_opt | 'UniStepXX' | 'VarStepXX' where XX is a number.
%              (defaut is UniStep1) E.g. UniStep5 will use piecewise
%              constant input functions with 5 different values  and a
%              constant time step depending on the simulation time.
%              input_opt is a structure with fields type ('UniStep' or
%              'VarStep') and cp (vector indicating number of control point
%              for each input signal)
%

% - pu0       default values for input parameters (default to 0)

%% default arguments
if ~exist('signals', 'var')||isempty('signals')
    signals='logged';
end

if ~exist('params', 'var')||isempty('params')
    params={};
end

if ~exist('p0', 'var')||isempty('p0')
    p0=[];
end

if ~exist('inputfn', 'var')||isempty('inputfn')
    inputfn='UniStep1';
end

if ~exist('pu0', 'var')||isempty('pu0')
    pu0=[];
end

if ~exist('options', 'var')||isempty('options')
    options=struct;
end


%% Copy the model into model_breach

% Give it a name
mdl_breach = [mdl '_breach'];
load_system(mdl);
close_system(mdl_breach,0);
save_system(mdl,mdl_breach);
close_system(mdl,0);
load_system(mdl_breach);

%% Checks options

% add scaling parameters to the lookup tables
if isfield(options, 'scale_lookup_tables')
    add_scale_param_lookup_tables(mdl_breach);
end

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

% Solver pane - times

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
cs.set_param('TimeSaveName', 'tout');   % Time
cs.set_param('SaveTime', 'on');   % Time
 
% Data Import/Export pane
cs.set_param('ExternalInput', '[t__, u__]');   % Input
cs.set_param('InspectSignalLogs', 'off');   % Record and inspect simulation output
cs.set_param('OutputSaveName', 'yout');   % Output
cs.set_param('ReturnWorkspaceOutputsName', 'out');   % Save simulation output as single object
cs.set_param('SaveCompleteFinalSimState', 'off');   % Save complete SimState in final state
cs.set_param('SaveFormat', 'StructureWithTime');   % Format
cs.set_param('SignalLoggingName', 'logsout');   % Signal logging name

if (~verLessThan('matlab','R2011a'))
    cs.set_param('DSMLoggingName', 'dsmout');   % Data stores logging name
    cs.set_param('SignalLoggingSaveFormat', 'Dataset');   % Signal logging format
    simfn = 'sim_breach';
else
    simfn = 'sim_breach_old_ver';   % call a different sim function for older versions of Simulink
                                    % FIXME: this likely does not work with 2010b...    
end

%% Find and Log input signals

lines = find_system(mdl_breach,'SearchDepth',1, 'FindAll', 'on', 'type', 'line');

for k = 1:numel(lines)
    bls = get(lines(k), 'SrcBlockHandle');
    if (bls ~= -1)
        type_bls = get(bls,'BlockType');
        
        % tests if we get an input
        if strcmp(type_bls, 'Inport')
            set(lines(k),'DataLoggingName', 'Use signal name', 'DataLogging',1 ,'Name',get(bls,'Name'));
        end
    end
end

%% Define Breach signals

% by default we log all input ports and all output ports

%% outputs
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

%% inputs
% find input ports and log their signals
sig_in ={};
InputOpt=[];
find_input_signals();

%% logged signals, by default none
sig_log = {};
find_logged_signals();


%% Scope signals
sig_scopes = find_scope_signals(mdl_breach);

%% define parameters
if exist('params','var')
    if (isempty(params))
        exclude = {'tspan','u__','t__'};
        assignin('base','tspan', 0:1);
        [params, p0] = filter_vars(mdl_breach, exclude);
    end
else
    exclude = {'tspan','u__','t__'};
    assignin('base','tspan', 0:1);
    [params, p0] = filter_vars(mdl_breach, exclude);
end

if isempty(p0)
    p0 = zeros(1,numel(params));
end

params = {params{:} U.params{:}};

if ~isempty(pu0)
    pu(1:numel(pu0)) = pu0;
end

%% Run the model for time 0 to check proper initialization and collect signal names
tspan = evalin('base', 'tspan;');
assignin('base','tspan',[0 0]);
assignin('base','t__',0);
assignin('base','u__',zeros(1, numel(sig_in)));

simout = sim(mdl_breach);
assignin('base','tspan',tspan);
[~,~, signals] = simout2X(simout);

%% Create the Breach structure
p0 = [zeros(1,numel(signals)) p0 pu];
Sys = CreateSystem(signals, params, p0'); % define signals and parameters

Sys.DimU = numel(sig_in);
Sys.InputList= sig_in;

if (exist('init_u'))
    Sys.init_u = init_u;
    Sys.InputOpt = InputOpt;
    [~, idx_inputs] = FindParamsInput(Sys);
    Sys.InputOpt.idx= idx_inputs;
end

Sys.type= 'Simulink';
eval(['Sys.sim = @' simfn ';']);
Sys.mdl= [mdl '_breach'];
Sys.Dir= pwd;
Sys.tspan = 0:t_step:t_end;
Sys.name = Sys.mdl;  % not great..

save_system(mdl_breach);
close_system(mdl_breach);

    function find_input_signals
        
        ins = find_system(mdl_breach,'SearchDepth',1, 'BlockType', 'Inport');
        
        for i = 1:numel(ins)
            
            in =  get_param(ins(i),'Name');
            in = regexprep(in,'\W','_');
            
            %Ensures port and its output line have the same name
            lh = get_param(ins(i), 'LineHandles');
            lh=lh{1}.Outport;
            set(lh, 'Name', in{1});
            sig_in = {sig_in{:} in{1}};
        end
        
        if (~isempty(sig_in))
            [c ia inew] = unique(sig_in);
            % reorder ...
            ia = sort(unique(ia));
            sig_in = sig_in(ia);
            
            % define the input generation function
            if ~exist('inputfn','var')
                % by default, constant 0
                Sys.InputOpt = struct('type','UniStep','cp', ones(1, sig_in));
                init_u =  @(sig_in,pts,tspan)UniStepSimulinkInput(1,sig_in,pts,tspan);
            else
                % backward compatibility with syntax UniStepXX and
                % VarStepXX (later might not work though)
                if ischar(inputfn)
                    init_u = [];
                    
                    pref = 'UniStep';
                    if regexp(inputfn, [pref '[0-9]+'])
                        cp = str2num(inputfn(numel(pref)+1:end));
                        inputfn = struct('type','UniStep','cp', cp*ones(1, numel(sig_in)));
                    else
                        pref = 'VarStep';
                        if regexp(inputfn, [pref '[0-9]+'])
                            cp = str2num(inputfn(numel(pref)+1:end));
                            inputfn = struct('type','VarStep','cp', cp*ones(1, numel(sig_in)));
                        end
                    end
                end
                if isstruct(inputfn)
                    
                    InputOpt = inputfn;
                    % checks number of control points
                    DimU = numel(sig_in);
                    if isscalar(InputOpt.cp)
                        InputOpt.cp = InputOpt.cp*ones(1, DimU);
                    end
                    
                    Sys.InputOpt = InputOpt;
                    Sys.InputList= sig_in;
                                       
                    U.params = FindParamsInput(Sys);
                    pu = zeros(1,numel(U.params));
                    init_u = @SimulinkInput;
                    
                    % checks interpolation method
                    
                    if ~isfield(InputOpt,'method')
                        InputOpt.method= 'previous';
                    end
                    
                    if ischar(InputOpt.method)
                        InputOpt.method = {InputOpt.method};
                    end
                    
                    
                    if numel(InputOpt.method)==1
                        method = InputOpt.method{1};
                        InputOpt.method = cell(1,DimU);
                        for iu = 1:DimU
                            InputOpt.method{iu} = method;
                        end
                    elseif numel(InputOpt.method)~= DimU
                        error('find_input_signal:invalid_interp_methods_number', ['Invalid number of interpolation methods, should be ' num2str(DimU)]);
                    end
                else
                    error('Invalid input function or input parameters.');
                end
            end
            
        else
            pu = [];
            U.params = {};
            cs.set_param('LoadExternalInput', 'off');   % Input
        end
        
        
    end

    function find_logged_signals
        
        if ~exist('signals','var')
            signals = {};
        end
        
        if ischar(signals)
            if strcmp(signals,'logged')
                try
                    [sig_log] = find_signals(mdl_breach);
                catch
                    error(['Could not analyze signals in ' mdl '. Check that your model is ' ...
                        'compiling and running.' ]);
                end
            end
            
        elseif iscell(signals)
            
            for k =  1:numel(signals);
                sig = signals{k};
                l = find_system(mdl_breach,'FindAll','on', 'type','line', 'Name', sig);
                if ~isempty(l)
                    %ensure that it is logged
                    set(lines(k),'DataLoggingName', 'Use signal name', 'DataLogging',1 ,'Name',sig);
                else
                    error(['Signal ' l 'not found in model.'] );
                end
            end
        end
    end

end
