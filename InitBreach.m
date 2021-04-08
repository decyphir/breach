function InitBreach(br_dir, force_init, varargin)
% InitBreach  initializes Breach, in particular adding paths to Breach directories


%% checks if global configuration variable is defined
global BreachGlobOpt

if nargin<2
    force_init = false;
end

if ~exist('br_dir', 'var')||isempty(br_dir)
    br_dir = which('InstallBreach');
    br_dir = fileparts(br_dir);    
end

if ~force_init && isfield(BreachGlobOpt, 'breach_dir')
    if  isequal(BreachGlobOpt.breach_dir, br_dir)
        return; % OK InitBreach has been run before
    end
end

addpath([br_dir filesep 'Core' filesep 'm_src']);

opt.verbose = 1;
opt.init_from_breachflows= false;
opt = varargin2struct_breach(opt, varargin{:});

addpath([br_dir filesep 'Core' filesep 'm_src']);
opt.verbose = 1;
opt.init_from_breachflows= false;
opt = varargin2struct_breach(opt, varargin{:});


%% remove old path, if any
br_all_dir = which('InstallBreach', '-all');
nb_dir = numel(br_all_dir);   
if  nb_dir>1
    for idir =  2: nb_dir
        br_old_dir =  fileparts(br_all_dir{idir});        
        all_paths = strsplit(path,pathsep);
        nb_paths = numel(all_paths);
        rm_path_list = {};
        for ii = 1:nb_paths
            if regexp(all_paths{ii},br_old_dir)
                %disp(['              ' all_paths{ii}]);
                rm_path_list = [rm_path_list all_paths{ii}];
            end
        end
        if ~isempty(rm_path_list)&&opt.verbose>=1
            disp(['Warning: removing paths in ' br_old_dir]);
            rmpath(rm_path_list{:});
            disp(' ');
        end
    end
end
    
%%  Make sure ModelsData exist
if ~exist( [br_dir filesep 'Ext' filesep 'ModelsData'], 'dir')
    mkdir([br_dir filesep 'Ext' filesep 'ModelsData']);
end


%%  Make sure ModelsData/ParallelTemp exist for parallel computing
if ~exist( [br_dir filesep 'Ext' filesep 'ModelsData' filesep 'ParallelTemp'], 'dir')
    mkdir([br_dir filesep 'Ext' filesep 'ModelsData' filesep 'ParallelTemp']);
end
%% Init
if opt.verbose>=1
    disp(['Initializing Breach from folder ' br_dir '...']);
end

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
    [br_dir filesep 'Core' filesep 'OutputGen'], ...
    [br_dir filesep 'Core' filesep 'ParamGen'], ...
    [br_dir filesep 'Core' filesep 'Diagnostics'], ...    
    [br_dir filesep 'Core' filesep 'Gui'], ...
    [br_dir filesep 'Core' filesep 'Gui' filesep 'elems'], ...
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
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'jsonlab'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'YAMLmatlab'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'YAMLmatlab' filesep 'external'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'YAMLmatlab' filesep 'extras'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'snobfit'], ...
    [br_dir filesep 'Ext' filesep 'Toolboxes' filesep 'minq'], ...
    };

addpath(list_path{:});
cd(cdr);

%% Lookfor extensions
if exist('InitBreachFlows', 'file')&&~opt.init_from_breachflows
    InitBreachFlows;
end


%% Init BreachGlobOpt global configuration variable

 % Some global constants
if ~isfield(BreachGlobOpt, 'MaxNumSamples')
        BreachGlobOpt.MaxNumSamples=100000;
end

% STL Initialization
BreachGlobOpt.disable_robust_linear_interpolation = 1;

if ~isfield(BreachGlobOpt,'RobustSemantics')
    BreachGlobOpt.RobustSemantics = 0;
end

if isfield(BreachGlobOpt, 'STLDB')
    if ~strcmp(class(BreachGlobOpt.STLDB), 'containers.Map')
        BreachGlobOpt.STLDB = containers.Map();
    end
else
    BreachGlobOpt.STLDB = containers.Map();
end

%% Store path list
BreachGlobOpt.breach_dir = br_dir;
BreachGlobOpt.list_path = list_path;

%% Force refreshing sl_customization menu % disabled fixed/proven useful
%sl_refresh_customizations;

warning('on',id);

end