%% Multi-objective and Optimization Constraints Support

%%
% Demonstrates the support for constraints in optimization 
% problems and multi-objective optimization. 

%% Introduction
% In Breach, there are two mechanisms to support constraints 
% and multi-objective optimization:
%
% * BreachRequirements object with multiple requirements 
% * BreachRequirements objects with pre-conditions
% 
% We first illustrate multi-objective using the automatic transmission
% system example. 
BrDemo.InitAutotrans

%% Multi-requirement 
% Multi-objective requirements can be created by creating a BreachRequirement object with more than one STL formula.  
req1 = STL_Formula('phi_speed', 'not (alw_[35,40] (speed[t]<speed_up and speed[t]>speed_low))'); 
req2 = STL_Formula('phi_rpm', 'not (alw_[35,40] (RPM[t]<rpm_up and RPM[t]>rpm_low))');

R = BreachRequirement({req1, req2}); 
R.SetParam({'speed_low','speed_up'}, [80, 85]);
R.SetParam({'rpm_low','rpm_up'}, [3000, 3100]);
R.PrintAll


%% Falisfication Problem Using Default Solver 
% The default solver has primitive multi-objective support. It has the same
% behavior as with one requirement being the conjunction of all
% requirements. However, we get display information for each of them. 
pb = FalsificationProblem(BrAutotrans_3ranges, R);
pb.solve();

%% Plotting 
% On plotting, the default view shows a bar plot with the number of
% satisfied and falsified requirements. Here the last sample falsified one
% requirement. Context menu allows to switch view to display robustness. 
Rlog = pb.GetLog();
BreachSamplesPlot(Rlog);

%% BreachRequirement with Pre-condition example
% Another mechanism to handle multiple requirements is to use
% pre-conditions within the BreachRequirement class. Pre-conditions are
% requirements that need to be satisfied to trigger the evaluation of the
% main requirements. We illustrate their use with the AFC model.
% 
BrDemo.InitAFC;
STL_ReadFile('AFC_simple_spec.stl');

%% 
% We define an input generator and 10 samples for the AFC model.  
B = BrAFC.copy();
InputGen = fixed_cp_signal_gen({'Pedal_Angle','Engine_Speed'}, [3 1], 'previous');
B.SetInputGen(InputGen);
B.SetParamRanges({'Pedal_Angle_u1', 'Pedal_Angle_u2', 'Engine_Speed_u0'}, [ 0 80; 0 80; 900 1100]);
B1 = B.copy();
B1.QuasiRandomSample(10);
B2 = B1.copy(); B3 = B.copy();


%% BreachRequirement Without Pre-condition
% We define a simple requirement evaluated on 10 traces.
R1 = BreachRequirement(AF_overshoot_req);
[v, V] = R1.Eval(B1)

%% BreachRequirement With Pre-conditions
% We now define a BreachRequirement object with pre-conditions on the
% input. 
R2 = BreachRequirement(AF_overshoot_req);
R2.SetParam('tol', 0.05);
R2.AddPreCond({'alw (Pedal_Angle[t]<60)',...
    'Pedal_Angle_u0 <= Pedal_Angle_u1',...
    'Pedal_Angle_u1 <= Pedal_Angle_u2'});
R3 = R2.copy();
B2.verbose = 1; % gives more information during simulation in particular 

%% 
% One of the motivation of pre-conditions is to make it possible to skip
% simulations when pre-conditions are input conditions and they are not
% satisfied. This happens when pre-conditions only involve input signals
% such as is the case here. 


%% Evaluation of a BreachRequirement With Pre-conditions
[v, V, Vp] = R2.Eval(B2) % third output are pre conditions satisfaction values



%% Solve problem with Input constraints 
% During falsification, when preconditions are not satisfied, the solver
% will use their satisfaction value as objective to maximize, actively
% trying to satisfy them. 
pb = FalsificationProblem(B3, R3); pb.freq_update = 5; pb.solve();

%% Plotting 
% When getting the logged samples, Breach indicates that traces were not
% computed because input signals did not satisfy the pre-condition(s). 
Rlog = pb.GetLog();
BreachSamplesPlot(Rlog);

