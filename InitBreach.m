function InitBreach(br_dir, force_init)
% InitBreach  initializes Breach, in particular adding paths to Breach directories

%% checks if global configuration variable is defined
global BreachGlobOpt

if ~exist('br_dir', 'var')||isempty(br_dir)
    br_dir = which('InstallBreach');
    br_dir = fileparts(br_dir);
end

if nargin<2
    force_init = false;
end

if ~force_init && isfield(BreachGlobOpt, 'breach_dir')
    if  isequal(BreachGlobOpt.breach_dir, br_dir)
        return; % OK InitBreach has been run before
    end
end

%% remove old path, if any
br_all_dir = which('InstallBreach', '-all');
nb_dir = numel(br_all_dir);   
if  nb_dir>1
    for idir =  2: nb_dir
        br_old_dir =   fileparts(br_all_dir{idir});
        all_paths = strsplit(path, ';');
        nb_paths = numel(all_paths);
        disp(['Warning: removing paths in ' br_old_dir]);
        rm_path_list = {};
        for ii = 1:nb_paths
            if strcmp(all_paths{ii}(1:min(numel(br_old_dir),end)),br_old_dir)
                %disp(['              ' all_paths{ii}]);
                rm_path_list = [rm_path_list all_paths{ii}];
            end
        end
        rmpath(rm_path_list{:});
        disp(' ');
    end
end
    
%%  Make sure ModelsData exist
if ~exist( [br_dir filesep 'Ext' filesep 'ModelsData'], 'dir')
    mkdir([br_dir filesep 'Ext' filesep 'ModelsData']);
end

%% Init
disp(['Initializing Breach from folder ' br_dir '...']);

id = 'MATLAB:dispatcher:nameConflict';
warning('off',id);

cdr = pwd;
cd(br_dir);

list_path = { ...
    br_dir, ...
    [br_dir filesep 'Core'], ...
    [br_dir filesep 'Core' filesep 'Init'], ...
    [br_dir filesep 'Core' filesep 'm_src'], ...
    [br_dir filesep 'Core' filesep 'Algos'], ...
    [br_dir filesep 'Core' filesep 'SignalGen'], ...
    [br_dir filesep 'Params'], ...
    [br_dir filesep 'Params' filesep 'm_src'], ...
    [br_dir filesep 'Params' filesep 'm_src' filesep 'sobolqr'], ...
    [br_dir filesep 'Params' filesep 'm_src' filesep 'niederreiter2'], ...
    [br_dir filesep 'Plots'], ...
    [br_dir filesep 'Plots' filesep 'm_src'], ...
    [br_dir filesep 'Online' filesep 'm_src'], ...
    [br_dir filesep 'Online' filesep 'bin'], ...
    [br_dir filesep 'Online' filesep 'simulink_stlib'], ...
    [br_dir filesep 'Examples'], ...
    [br_dir filesep 'Ext' filesep 'Models'], ...
    [br_dir filesep 'Ext' filesep 'ModelsData'], ...
    [br_dir filesep 'Ext' filesep 'Specs'], ...
    [br_dir filesep 'Ext' filesep 'Classes'], ...
    [br_dir filesep 'Ext' filesep 'Specs' filesep 'STLib'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'optimize'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'DataHash'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'simulink_custom'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'sundials' filesep 'sundialsTB' ], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'sundials' filesep 'sundialsTB' filesep 'cvodes'], ...
    };

addpath(list_path{:});

%% Init BreachGlobOpt options and fourre-tout global variable
if exist('BreachGlobOpt.mat')
    load BreachGlobOpt;
    
    % Convert BreachGlobOpt into global
    BreachGlobOptTmp = BreachGlobOpt;
    clear BreachGlobOpt;
    global BreachGlobOpt;
    BreachGlobOpt = BreachGlobOptTmp;
    clear BreachGlobOptTmp;
    BreachGlobOpt.RobustSemantics = 0 ; % 0 by default, -1 is for left time robustness, +1 for right, inf for sum ?
    
else
    
    BreachGlobOpt.breach_dir = br_dir;
    
    if ~isfield(BreachGlobOpt,'RobustSemantics')
        BreachGlobOpt.RobustSemantics = 0;
    end
    
end
cd(cdr);


%% Init STL_Formula database

if isfield(BreachGlobOpt, 'STLDB')
    if ~strcmp(class(BreachGlobOpt.STLDB), 'containers.Map')
        BreachGlobOpt.STLDB = containers.Map();
    end
else
    BreachGlobOpt.STLDB = containers.Map();
end

%% Store path_list for when we want to remove it
BreachGlobOpt.list_path = list_path;

%% Force refreshing sl_customization menu
sl_refresh_customizations;

warning('on',id);

end