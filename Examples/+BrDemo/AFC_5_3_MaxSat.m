%% Maximizing Satisfaction Problems

%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
BrDemo.InitAFC;
BrAFC

%%
% Load some properties.
tol = 6.6e-3;
R = BreachRequirement('AFC_simple_spec.stl', 'AF_alw_ok');
R.SetParam('tol', tol);


%% Maximizing Satisfaction
% Instead of minimizing satisfaction for falsification, we might want to
% maximize it to make it positive and as robust as possible.

%% 
% One possible use case is tuning control parameters.
B = BrAFC.copy();
% Create the max sat problem and solve it
max_sat_pb = MaxSatProblem(B, R,...
                                {'ki', 'kp'}, ...       % we look for a better PI controller. 
                                [0. 0.1; 0.01 0.1]);
max_sat_pb.solve();

%%
% The best robustness value obtained is positive, which means we found
% 'better' parameters for PI controller that satisfy the specification.

%% Maximizing Satisfaction - Plot
% Maybe it is not so good though.

AFC_Best = max_sat_pb.GetBrSet_Best();
AFC_Best.BrSet.PlotSignals({'AF'}, [], {'LineWidth', 2});
set(gca, 'XLim', [10 40], 'FontSize',14, 'LineWidth', 2);
plot([0 41], (1+tol)*[14.7 14.7],'r');
plot([0 41], (1-tol)*[14.7 14.7],'r');

%% Adjusting Requirement
% We maximized the satisfaction of AF_alw_tol, but this resulted in an
% undesired oscillatory behavior. To fix this, we use a new formula that
% requires that AF settles for at least 1s every 10s.
STL_ReadFile('AFC_settling_spec.stl')
type AFC_settling_spec.stl

%% Maximizing for the New Requirement
% Defining and solving the max sat problem with settling property

B = BrAFC.copy();
R = BreachRequirement('AF_alw_settle and AF_alw_ok');
R.SetParam('tol', tol);
R.SetParam('dt', 0.1);
R.SetParam('epsi', 1e-2);
R.SetParam('t_start', 10);
R.SetParam('t_end', 10);

% Create the max sat problem and solve it
max_sat_pb = MaxSatProblem(B, R,...
                           {'ki', 'kp'}, ...       
                           [0. 0.1; 0.01 0.1]);
max_sat_pb.solve();

%%
% The solver found a satisfactory solution that should behave better. 

%% Maximizing for the New Requirement - Plot

AFC_Best = max_sat_pb.GetBrSet_Best();
AFC_Best.BrSet.PlotSignals({'AF'}, [], {'LineWidth', 2});
set(gca, 'XLim', [10 40], 'FontSize',14, 'LineWidth', 2);
plot([0 41], (1+tol)*[14.7 14.7],'r');
plot([0 41], (1-tol)*[14.7 14.7],'r');

%% Maximizing Satisfaction For Multiple Inputs
% Tuning PI for only one input is likely not sufficient. Below we try to
% find a PI controller for a set of 9 inputs.

AFC_Grid = BrAFC.copy();
AFC_Grid.SetParamRanges({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [10, 20; 10 60]);
AFC_Grid.GridSample([3 3]);

max_sat_problem_grid = MaxSatProblem(AFC_Grid, AF_alw_settle, {'ki', 'kp'}, [0. 0.2; 0.0 0.05]);                               
max_sat_problem_grid.solve();

%% 
% If the solver doesn't find a solution, we can again try another one,
% increase the simulation budget, change the problem, etc. 

AFC_Best = max_sat_problem_grid.GetBrSet_Best();
AFC_Best.BrSet.PlotSignals({'AF'}, [], {'LineWidth', 2});
set(gca, 'XLim', [10 40], 'FontSize',14, 'LineWidth', 2);
plot([0 41], (1+tol)*[14.7 14.7],'r');
plot([0 41], (1-tol)*[14.7 14.7],'r');
