function [BrSynth, SynthProb, FalsifProb] = BreachCEGIS(SynthProb, FalsifProb, options_in)
%CEGIS Implements a Counter-Example Guided Inductive Synthesis strategy
%


%% processing options
% if zero argument, returns an option structure with defaults
if nargin==0
    BrSynth= struct('iter_max', 10);
    return
end

% no option, use defaults
if nargin == 2
    options_in = struct();
end

% option provided, make sure all fields are initialized 
options = BreachCEGIS();
opt_in_fields = fieldnames(options_in);
for  ifld=1:numel(opt_in_fields)
    options.(opt_in_fields{ifld}) = options_in.(opt_in_fields{ifld});
end

% assigning option parameters
iter_max = options.iter_max;

%% Main loop
iter = 1;
cont = true;
while (cont)
    clc;
    %% Synthesis step
    fprintf('Iter %d/%d\n', iter,iter_max)
    fprintf('Synthesis step\n');
    fprintf('--------------\n');
    SynthProb.solve();
    BrSynth = SynthProb.GetBrSet_Best();
    
    %% Falsification step
    fprintf('Counter-Example step\n');
    fprintf('------------------------\n');
    FalsifProb.BrSet = BrSynth;
    FalsifProb.ResetObjective();
    FalsifProb.solve();
    BrFalse = FalsifProb.BrSet_False;
    if isempty(BrFalse)
        return
    end
    
    %% Update parameter synthesis problem
    SynthProb.BrSet.Concat(BrFalse);
    SynthProb.ResetObjective();
    SynthProb.BrSys.Sys.Verbose=0;
    
    iter = iter+1;
    cont = iter<iter_max;
end


end