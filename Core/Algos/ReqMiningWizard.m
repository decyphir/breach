function rq_pb= ReqMiningWizard(varargin)

options.ReqMiningProblem = []; 
options.ParamSet = [];
options.Requirement = [];
opt = varargin2struct_breach(options, varargin{:});

%% Step 1: Choose existing problem or create new one from set and requirement
% Detect existing sets and problems

% Falsif problem
reqmining_pbs = get_vars_in_base('ReqMiningProblem');
if isempty(opt.ReqMiningProblem)
    opt.ReqMiningProblem = 'New Problem';
end

choices.ReqMiningProblem = ['New Problem' reqmining_pbs];
tips.ReqMiningProblem = 'Choose to create a new requirement mining problem or pick an existing one.';

% Parameter set 
breach_sets = get_vars_in_base('BreachSet');
if isempty(opt.ParamSet)
    opt.ParamSet = breach_sets{1};
end
breach_sets = [{''} breach_sets];
choices.ParamSet = breach_sets;
tips.ParamSet = 'Pick a parameter set. If empty, reuse from existing problem. The parameter set should have both system and requirement variables.';

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
tips.Requirement = 'Choose a requirement for falsification and mining. If isempty, reuse from existing problem.';

% Solver choice
opt.FalsifSolver = 'global_nelder_mead';
choices.FalsifSolver = {'basic', 'global_nelder_mead', 'simulannealbnd', 'cmaes'};
tips.FalsifSolver = 'Choose a solver for the falsification steps.';

% Number of iterations
opt.NumFalsifSteps = 10;
choices.NumFalsifSteps = 'int';
tips.NumFalsifSteps = 'Maximum number of falsification steps after requirement is updated.';

% Checks if there is something to chew on
if isempty(reqmining_pbs)&&(isempty(breach_sets)||isempty(reqs))
    error('Load or create a parameter set and a requirement first.');
end

gu1 = BreachOptionGui( 'ReqMining Problem', opt, choices, tips);
uiwait(gu1.dlg);

if isempty(gu1.output)
    rq_pb = [];
    return;
end

% Create falsif problem if needed
if isequal(gu1.output.ReqMiningProblem,'New Problem');
    B = evalin('base', gu1.output.ParamSet);
    if B.Specs.isKey(gu1.output.Requirement)
        req = B.Specs(gu1.output.Requirement);
    else
        req = evalin('base', gu1.output.Requirement);
    end
    rq_pb = ReqMiningProblem(B , req );
else
    rq_pb = evalin('base', gu1.output.ReqMiningProblem);
end

rq_pb.iter_max = gu1.output.NumFalsifSteps; 

%% Step 2 Solver options

opt = struct; choices = struct; tips = struct;
rq_pb_name = rq_pb.whoamI; 
if isequal(rq_pb_name, '__Nobody__')
    rq_pb_name = evalin('base',   'matlab.lang.makeUniqueStrings(''rq_pb'', who);');
end


set_opt('ProblemName',...
            rq_pb_name,...
            'string',...
            'Name for a variable the workspace storing the ReqMiningProblem object.'); 

set_opt('use_parallel',...
            rq_pb.falsif_pb.use_parallel,...
            'bool',...
            'If true, the solver will be allowed to use parallel computation if it supports it.'); 

set_opt('max_time', ...
    rq_pb.falsif_pb.max_time,...
    'int',...
    'Pick a total computational time (seconds) after which the solver will stop looking for a solution.');         
        
set_opt('max_obj_eval', ...
    rq_pb.falsif_pb.max_obj_eval,...
    'int',...
    'Pick a total objective function evaluation number, i.e., number of simulations, after which the solver will stop looking for a solution.');         

set_opt('verbose',...
    rq_pb.falsif_pb.verbose,...
    'int',...
    'Enter a number greater or equal to 0. The larger the noiser.');

gu2 = BreachOptionGui( 'ReqMining Problem', opt, choices, tips);

%% Step 3 choose solver specific options 
[solver_opt, is_gui] =  rq_pb.falsif_pb.setup_solver(gu1.output.FalsifSolver, true);
%fp.solve();
if is_gui
    gu2.merge_options(solver_opt.output, solver_opt.choices, solver_opt.tips);
end

uiwait(gu2.dlg);
if isempty(gu2.output)
    rq_pb = [];
    return;
end
rq_pb.falsif_pb.use_parallel = gu2.output.use_parallel; 
rq_pb.falsif_pb.max_time = gu2.output.max_time;
rq_pb.falsif_pb.max_obj_eval = gu2.output.max_obj_eval;
rq_pb.falsif_pb.verbose = gu2.output.verbose;

% solver option
if is_gui
    rq_pb.falsif_pb.solver_options.use_param_set_as_init = gu2.output.use_param_set_as_init;
    rq_pb.falsif_pb.solver_options.start_at_trial = gu2.output.start_at_trial;
    rq_pb.falsif_pb.solver_options.nb_new_trials = gu2.output.nb_new_trials;
    rq_pb.falsif_pb.solver_options.nb_local_iter = gu2.output.nb_local_iter;
end

assignin('base', gu2.output.ProblemName, rq_pb);


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
