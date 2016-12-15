InitBreach;

max_rob = 10000;              % initial robustness is between -max_rob, +max_rob
rob_up_lim = -max_rob;  % stops when rob_up is below 
rob_low_lim = max_rob; 
rob_diff_lim = -1;

phi_autotrans = 'gear[t]>0';

% Model
mdl = 'Autotrans_online';
Br = BreachSimulinkSystem(mdl); 

% define the formula
phi1 = STL_Formula('phi1', '(alw (speed[t]<vmax)) and (alw (RPM[t]<rpm_max))');
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