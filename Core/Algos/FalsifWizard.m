function fp= FalsifWizard(varargin)

options.FalsificationProblem = []; 
options.ParamSet = [];
options.Requirement = [];
opt = varargin2struct_breach(options, varargin{:});

%% Step 1: Choose existing problem or create new one from set and requirement
% Detect existing sets and problems

% Falsif problem
falsif_pb = get_vars_in_base('FalsificationProblem');
if isempty(opt.FalsificationProblem)
    opt.FalsificationProblem = 'New Problem';
end

choices.FalsificationProblem = ['New Problem' falsif_pb];
tips.FalsificationProblem = 'Choose to create a new falsification problem or pick an existing one.';

% Parameter set 
breach_sets = get_vars_in_base('BreachSet');
if isempty(opt.ParamSet)
    opt.ParamSet = breach_sets{1};
end
breach_sets = [{''} breach_sets];
choices.ParamSet = breach_sets;
tips.ParamSet = 'Pick a parameter set. If empty, reuse from existing problem.';

% Requirement
reqs = get_vars_in_base('STL_Formula');
reqs = [{''} reqs];
if ~isempty(opt.ParamSet)  % if a paramset is specified, make sure its specs are listed 
    try
        Breqs = evalin('base', [ opt.ParamSet '.Specs' ]);
        reqs = union(reqs, Breqs.keys);
    end
 end

if isempty(opt.Requirement)
    opt.Requirement = reqs{1};
end
assert(isempty(reqs{1}));
choices.Requirement = reqs;
tips.Requirement = 'Choose a requirement to falsify. If isempty, reuse from existing problem.';

% Solver choice
opt.Solver = 'global_nelder_mead';
choices.Solver = {'basic', 'global_nelder_mead', 'simulannealbnd', 'cmaes'};
tips.Solver = 'Choose a solver.';

% Checks if there is something to chew on
if isempty(falsif_pb)&&(isempty(breach_sets)||isempty(reqs))
    error('Load or create a parameter set and a requirement first.');
end

gu1 = BreachOptionGui( 'Falsification Problem', opt, choices, tips);
uiwait(gu1.dlg);

if isempty(gu1.output)
    fp = [];
    return;
end

% Create falsif problem if needed
if isequal(gu1.output.FalsificationProblem,'New Problem');
    B = evalin('base', gu1.output.ParamSet);
    if B.Specs.isKey(gu1.output.Requirement)
        req = B.Specs(gu1.output.Requirement);
    else
        req = evalin('base', gu1.output.Requirement);
    end
    fp = FalsificationProblem(B , req );
else
    fp = evalin('base', gu1.output.FalsificationProblem);
end

%% Step 2 Solver options

opt = struct; choices = struct; tips = struct;
fp_name = fp.whoamI; 
if isequal(fp_name, '__Nobody__')
    fp_name = evalin('base',   'matlab.lang.makeUniqueStrings(''falsif_pb'', who);');
end


set_opt('ProblemName',...
            fp_name,...
            'string',...
            'Name for a variable the workspace storing the FalsificationProblem object.'); 

set_opt('use_parallel',...
            fp.use_parallel,...
            'bool',...
            'If true, the solver will be allowed to use parallel computation if it supports it.'); 

set_opt('max_time', ...
    fp.max_time,...
    'int',...
    'Pick a total computational time (seconds) after which the solver will stop looking for a solution.');         
        
set_opt('max_obj_eval', ...
    fp.max_obj_eval,...
    'int',...
    'Pick a total objective function evaluation number, i.e., number of simulations, after which the solver will stop looking for a solution.');         

set_opt('verbose',...
    fp.verbose,...
    'int',...
    'Enter a number greater or equal to 0. The larger the noiser.');

gu2 = BreachOptionGui( 'Falsification Problem', opt, choices, tips);

%% Step 3 choose solver specific options 
[solver_opt, is_gui] =  fp.setup_solver(gu1.output.Solver, true);
%fp.solve();
if is_gui
    gu2.merge_options(solver_opt.output, solver_opt.choices, solver_opt.tips);
end

uiwait(gu2.dlg);
if isempty(gu2.output)
    fp = [];
    return;
end
fp.use_parallel = gu2.output.use_parallel; 
fp.max_time = gu2.output.max_time;
fp.max_obj_eval = gu2.output.max_obj_eval;
fp.verbose = gu2.output.verbose;

% solver option
if is_gui
    fp.solver_options.use_param_set_as_init = gu2.output.use_param_set_as_init;
    fp.solver_options.start_at_trial = gu2.output.start_at_trial;
    fp.solver_options.nb_new_trials = gu2.output.nb_new_trials;
    fp.solver_options.nb_local_iter = gu2.output.nb_local_iter;
end

assignin('base', gu2.output.ProblemName, fp);


%% Step 4 

    function set_opt(opt_name, opt_value, opt_choices, opt_tip)
        opt.(opt_name) = opt_value;
        choices.(opt_name) = opt_choices;
        tips.(opt_name) = opt_tip;
    end

end

function sets = get_vars_in_base(type)
sets = {};
ws_var = evalin('base', 'who');
for iv= 1:numel(ws_var)
    % is this a BreachSet?
    varname = ws_var{iv};
    if  ~isequal(varname, 'ans')
        BB__ = evalin('base', varname );
        if isa(BB__, type) % found one, keep it
            sets  = [sets  varname];
        end
    end
end
end
