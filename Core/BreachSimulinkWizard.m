function Br = BreachSimulinkWizard(mdl, Br, varargin)
% Takes a model mdl, and through a sequence of GUIs, creates an interface.

if ~exist('Br', 'var')
    Br = [];
end


%% Arguments parsing and sanity check
options.RunFromSimulink = 0;
try 
    options.VariableName =Br.whoami;
catch
    options.VariableName = 'Br';
end
options.FindTables = false;
options.MaxNumTabParam = 10;
options.Verbose = 1;

if ~isempty(Br)
    options.InitFn = Br.InitFn;
else
    options.InitFn = '';
end

if nargin >=3
    options = varargin2struct_breach(options, varargin{:});
end

if nargin==0
    try
        fprintf_log('Finding current simulink model...')
        mdl = bdroot;
        fprintf_log('current model is %s\n', mdl);
    catch
        fprintf_log('failed.\n')
        error('BreachSimulinkWizard:no_open_model_not_found',[mdl ' does not reference a Simulink model']);
    end
else
    load_system(mdl);
end

if regexp(mdl, '.*_breach')
    y= input(['The current model appears to be a breach copy of  ' mdl(1:end-7) '.\nDo you want to proceed with interfacing ' mdl(1:end-7) '?']);
    if ~y
        Br = [];
        return;
    else
        mdl = mdl(1:end-7);
        load_system(mdl);
    end
end
%% user options

opt = rmfield(options, 'RunFromSimulink');

choices.VariableName = 'string';
tip.VariableName = 'Enter a name for the BreachSimulinkSystem object variable.'; 
choices.Verbose = 'int';
tip.Verbose = 'Enter a number greater or equal to 0. The larger the noiser.';
choices.FindTables = 'bool';
tip.FindTables = 'If true, Breach will create a parameter for each table elements.';
choices.MaxNumTabParam = 'int';
tip.MaxNumTabParam = 'Enter a number greater or equal to 1. Breach will only create parameters for as many elements in a table (avoids creating too many parameters in case of large tables).';
choices.InitFn = 'string';
tip.InitFn = 'Enter the name of an initialization function, if one is needed. Leave empty otherwise. The function must be a script or a function without arguments executed in the base workspace.';

gu = BreachOptionGui(['Choose options for interfacing model ' mdl], opt,choices, tip);
uiwait(gu.dlg);

opt = gu.output;

if isempty(opt)
    Br=[];
    clean();
    return
end

assignin('base', 'opt__', rmfield(opt, 'VariableName'));
opt.InitFn = [];

%% Run initialization function if there is one
if ~isempty(opt.InitFn)
fprintf_log('Running initialization function  %s ... ',opt.InitFn);
try
    assignin('base', 'InitFn__', opt.InitFn);
    evalin('base',opt.InitFn);
catch
     fprintf_log('failed with message: %s.\n', err_msg);
     return;
end
end


%% Save model

fprintf_log('Saving %s ... ', mdl);
try
    save_system(mdl);
    fprintf_log('done.\n');
catch
    err_msg = lasterr;
    fprintf_log('failed with message: %s.\n', err_msg);
    return
end


%%  Execute model for 0 time (forces compilation, checks and stuff)
sim_cmd = ['sim(''' mdl ''',[0 0]);'] ;
fprintf_log('Trying command: ''%s'' ... ', sim_cmd);
crd = pwd;
cd(BreachGetModelsDataPath());
try  
    evalin('base', sim_cmd);
    fprintf_log('done.\n');
catch
    fprintf_log('failed.\n')
    help_msg = sprintf(['Error: To create an interface, ' ...
        'Breach requires that Simulink sim command runs ' ...
        'without error for model %s. Try running the model '...
        'in Simulink editor and fix any error.\n'  ...
        'Simulink returned the following:\n\n'], mdl);
    [err_msg] = lasterr;
    help_msg = sprintf([help_msg '%s'],  err_msg);
    title = 'Error with model compilation';
    disp_log(help_msg, title);
    clean();
    if options.RunFromSimulink
        open_system(mdl);
    end
    cd(crd);
    return;
end
cd(crd);

%%  Creates default interface
Name = get_param(mdl, 'Name');
close_system(mdl);
Ball = evalin('base', ['BreachSimulinkSystem(''' Name ''', ''all'', [], {}, [], opt__)']);
if isempty(Br)
    Br = Ball;
    assignin('base','IG__', []);
else
    assignin('base', 'IG__', Br.InputGenerator);
end

%%
signals_all = Ball.GetSignalList();
signals = Br.GetSignalList();
niou_sigs = select_cell_gui(signals_all, signals, 'Select signals from list');
if isequal(niou_sigs,0) % pushed cancel
    if options.RunFromSimulink % get back to Simulink editor
        open_system(mdl);
    end
    return;
end

assignin('base','sigs__', niou_sigs);
params_all = Ball.GetPlantParamList();
params = Br.GetPlantParamList();

%% Params in workspace
if ~isempty(params_all)
    params= select_cell_gui(params_all, params, 'Select parameters from list');
   
    if isequal(params,0) % pushed cancel
        if options.RunFromSimulink % get back to Simulink editor
            open_system(mdl);
        end
        return;
    end
end
 assignin('base', 'params__', params)
evalin('base', [opt.VariableName '=BreachSimulinkSystem(''' Name ''', params__,[],sigs__,IG__,opt__ );']);
Br = evalin('base',opt.VariableName);
gu = Br.SetInputGenGUI();
uiwait(gu);

clean;

if options.RunFromSimulink
    Br.RunGUI();
end

%% Auxiliary  functions


    function disp_log(st, title)
        % disp in command line or in a dialog box
        if options.RunFromSimulink
            g =  msgbox(st, title);
            uiwait(g);
        else
            disp(st);
        end
    end


%% log version of fprintf
    function fprintf_log(varargin)
        if options.Verbose
            fprintf(varargin{:});
        end
    end

%% cleaning
    function clean()
        evalin('base', 'clear params__ sigs__ opt__ IG__ InitFn__;');
    end

end