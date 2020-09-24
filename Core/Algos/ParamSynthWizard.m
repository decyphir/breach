function pb= ParamSynthWizard(varargin)

options.ParamSynthProblem = []; 
options.ParamSet = [];
options.Requirement = [];
opt = varargin2struct_breach(options, varargin{:});

%% Step 1: Choose existing problem or create new one from set and requirement
% Detect existing sets and problems

% Parameter Synthesis problem
param_synth_pb = get_vars_in_base('ParamSynthProblem');
if isempty(opt.ParamSynthProblem)
    opt.ParamSynthProblem = 'New Problem';
end

choices.ParamSynthProblem = ['New Problem' param_synth_pb];
tips.ParamSynthProblem = 'Choose to create a new falsification problem or pick an existing one.';

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
    Breqs = evalin('base', [ opt.ParamSet '.Specs' ]);
    reqs = union(reqs, Breqs.keys);
end

if isempty(opt.Requirement)
    opt.Requirement = reqs{1};
end
assert(isempty(reqs{1}));
choices.Requirement = reqs;
tips.Requirement = 'Choose a requirement for parameter synthesis. If isempty, reuse from existing problem.';


% Checks if there is something to chew on
if isempty(param_synth_pb)&&(isempty(breach_sets)||isempty(reqs))
    error('Load or create a parameter set and a requirement first.');
end

gu1 = BreachOptionGui( 'Parameter Synthesis Problem', opt, choices, tips);
uiwait(gu1.dlg);

if isempty(gu1.output)
    pb = [];
    return;
end

% Create falsif problem if needed
if isequal(gu1.output.ParamSynthProblem,'New Problem');
    B = evalin('base', gu1.output.ParamSet);
    if B.Specs.isKey(gu1.output.Requirement)
        req = B.Specs(gu1.output.Requirement);
    else
        req = evalin('base', gu1.output.Requirement);
    end
    pb = ParamSynthProblem(B , req );
else
    pb = evalin('base', gu1.output.ParamSynthProblem);
end

%% Step 2 Solver options

opt = struct; choices = struct; tips = struct;
fp_name = pb.whoamI; 
if isequal(fp_name, '__Nobody__')
    fp_name = evalin('base',   'matlab.lang.makeUniqueStrings(''param_synth_pb'', who);');
end


set_opt('ProblemName',...
            fp_name,...
            'string',...
            'Name for a variable the workspace storing the ParamSynthProblem object.'); 

set_opt('max_time', ...
    pb.max_time,...
    'int',...
    'Pick a total computational time (seconds) after which the solver will stop looking for a solution.');         
        
set_opt('max_obj_eval', ...
    pb.max_obj_eval,...
    'int',...
    'Pick a total objective function evaluation number, i.e., number of simulations, after which the solver will stop looking for a solution.');         

set_opt('verbose',...
    pb.verbose,...
    'int',...
    'Enter a number greater or equal to 0. The larger the noiser.');

gu2 = BreachOptionGui( 'ParamSynth Problem', opt, choices, tips);
uiwait(gu2.dlg);
pb.max_time = gu2.output.max_time;
pb.max_obj_eval = gu2.output.max_obj_eval;
pb.verbose = gu2.output.verbose;

assignin('base', gu2.output.ProblemName, pb);

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
