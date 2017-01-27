InitBreach;

max_rob = 10000;              % initial robustness is between -max_rob, +max_rob
rob_up_lim = 0;  % stops when rob_up is below 
rob_low_lim = 0; 
rob_diff_lim = 0;

phi_autotrans = '(alw_[0, 40] (speed[t]<160)) and (alw_[0,40] (RPM[t]<5000))';

% Model
mdl = 'Autotrans_online';
Br = BreachSimulinkSystem(mdl); 

% define the formula
phi1 = STL_Formula('phi1', 'alw (rob_up[t]>0)');
phi1 = set_params(phi1,{'vmax', 'rpm_max'}, [160 4500]);

% define the input parameters and ranges 
input_params.names = {'Throttle_u0'};
input_params.ranges = [0 100];
 
% defines the falsification problem and solve it
falsif_pb = FalsificationProblem(Br, phi1, input_params.names, input_params.ranges);

% solves using the default solver
falsif_pb.solve();

% collect the falsifying trace 
BrFalse = falsif_pb.GetBrSet_False() 
BrFalse.PlotRobustSat(phi1)