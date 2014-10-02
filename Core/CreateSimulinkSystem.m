function Sys = CreateSimulinkSystem(mdl, signals, params, p0,  inputfn, pu0 )
% CreateSimulinkSystem  Create a Breach system structure from a simulink model
%  
% Synopsis: Sys =  CreateSimulinkSystem(mdl, signals, params, p0, inputfn, pu0)
%
%  - mdl       name of the Simulink model - note that breach will use its own copy of this model
%  - signals   can be either 'logged' or {'s1','s2', ..., 'sn'} where s1, s2,..., sn are logged signals in the model  
%              if signals=='logged', Breach will find and use all logged signals in the model 
%  - params    {'p1',...,'pn'} where pi are tunable parameters in mdl. If
%              empty or absent, Breach will try to find tunable parameters
%  - p0        default values for the tunable parameters (0 if empty)
%  - inputfn   'UniStepXX' | 'VarStepXX' | 'UniPWAXX' where XX is a number (defaut is UniStep1)
%               E.g. UniStep5 will use piecewise constant input functions with 5 different values 
%               and a constant time step depending on the simulation time.  
%  - pu0        default values for input parameters (default to 0)

%% default arguments

switch nargin
    case 1
        signals ='logged';
        params={};
        p0=[];
        inputfn = 'UniStep1';
    case 2
        params={};
        p0=[];
        inputfn = 'UniStep1';
    case 3
        p0=[];
        inputfn = 'UniStep1';
    case 4    
        inputfn = 'UniStep1';
end


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
% cs.set_param('SaveOutput', 'on');   % Output 
  cs.set_param('SaveTime', 'on');   % Time 

% cs.set_param('SolverType', 'Fixed-step');   % Type 

% Solver pane
%  cs.set_param('FixedStep', 'auto');   % Fixed-step size (fundamental sample time) 
%  cs.set_param('Solver', 'ode5');   % Solver 
  cs.set_param('StartTime', '0.0');   % Start time 
  cs.set_param('StopTime', 'tspan(end)');   % Stop time 

% Data Import/Export pane
  cs.set_param('ExternalInput', '[t__, u__]');   % Input 
% cs.set_param('FinalStateName', 'xFinal');   % Final states 
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
  end
  cs.set_param('TimeSaveName', 'tout');   % Time 
    
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

  % outputs
  o = find_system(mdl_breach,'SearchDepth',1, 'BlockType', 'Outport');
 
  sig_out= {};
  for i = 1:numel(o)
    nm = regexprep(o(i),[mdl_breach '/'],'');    
    sig_out = {sig_out{:}, nm{:}};
  end
 
%% inputs 
  
  % find input ports and log their signals 
  
  ins = find_system(mdl_breach,'SearchDepth',1, 'BlockType', 'Inport');
  
  sig_in ={};
  for i = 1:numel(ins)
    in =  get_param(ins(i),'Name'); 
    sig_in = {sig_in{:} in{1}};
  end  
  
  if (~isempty(sig_in))
    [c ia inew] = unique(sig_in);
    % reorder ...
    ia = sort(unique(ia));
    sig_in = sig_in(ia);
     
    % define the input generation function
    
    if ~exist('inputfn')
      % by default, constant 0
      init_u =  @(sig_in,pts,tspan)UniStepSimulinkInput(1,sig_in,pts,tspan);
    else 
      if isstr(inputfn)
        
        init_u = [];
                
        pref = 'UniStep';  % FIXME find a way to set interpolate data for each input block      
        if regexp(inputfn, [pref '[0-9]+'])
          cp = str2num(inputfn(numel(pref)+1:end));
          init_u =  @(sig_in,pts,tspan)UniStepSimulinkInput(cp ,sig_in,pts,tspan);        
        end
      
        if isempty(init_u)
          pref = 'VarStep';        
          if regexp(inputfn, [pref '[0-9]+'])
            cp = str2num(inputfn(numel(pref)+1:end));
            init_u =  @(sig_in,pts,tspan)VarStepSimulinkInput(cp ,sig_in,pts,tspan);        
          end
        end
    
        if isempty(init_u) 
          pref = 'UniPWA';  % FIXME find a way to set interpolate data for each input block       
          if regexp(inputfn, [pref '[0-9]+'])
            cp = str2num(inputfn(numel(pref)+1:end));
            init_u =  @(sig_in,pts,tspan)UniPWASimulinkInput(cp ,sig_in,pts,tspan);        
          end
        end
                
        
        if isempty(init_u)
          eval(['init_u=@' inputfn ';']);
        end              
      else
        error('Invalid input function.');
      end
    end  
    
    U = init_u(sig_in,[],[]);
    pu  = U.p0;
    
  else
    pu = [];
    U.params = {};
    cs.set_param('LoadExternalInput', 'off');   % Input  
    
  end
  
  %% logged signals, by default none
  if ~exist('signals')
    signals = {};
  end

  sig_log = {};
  if isstr(signals)
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
  
  % finds duplicates : signals logged that are inputs or outputs
  keep = ones(1,numel(sig_log));
  for k = 1:numel(sig_log);
    for ks= 1:numel(sig_in)
      if strcmp(sig_log{k}, sig_in{ks}) 
        keep(k)= 0;
      end
    end
    
    for ks = 1:numel(sig_out)
      if strcmp(sig_log{k}, sig_out{ks}) 
        keep(k)= 0;
      end
    end
  end
  sig_log= sig_log(find(keep));
  
  signals = {sig_out{:} sig_log{:} sig_in{:}};
 
  %% TODO For all signals, check wether it's multi-dimensional
       
  
%% define parameters 
  
  if exist('params')
    if (isempty(params))  
      exclude = {'tspan','u__','t__'};
      assignin('base','tspan', 0:1);
      [params p0] = filter_vars(mdl_breach, exclude);      
    end
    
  else
     exclude = {'tspan','u__','t__'};
     assignin('base','tspan', 0:1);
     [params p0] = filter_vars(mdl_breach, exclude);  
  end   
  
  if ~exist('p0')
      p0 = zeros(1,numel(params));
  end
  
  if isempty(p0)
      p0 = zeros(1,numel(params));
  end
    
  params = {params{:} U.params{:}};
  
  if (exist('pu0'))  
      pu(1:numel(pu0)) = pu0;
  end
      
  p0 = [zeros(1,numel(signals)) p0 pu];
     
  Sys = CreateSystem(signals, params, p0'); % define signals and parameters

  Sys.DimY = numel(sig_out);
  Sys.DimU = numel(sig_in);  
  if (exist('init_u'))
    Sys.init_u = init_u;
  end
  
  Sys.type= 'Simulink';
  eval(['Sys.sim = @' simfn ';']); 
  Sys.mdl= [mdl '_breach'];
  Sys.Dir= pwd;
  
  save_system(mdl_breach);  
  close_system(mdl_breach);

  