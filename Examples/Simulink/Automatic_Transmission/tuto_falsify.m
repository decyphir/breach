InitBreach;

% Model
mdl = 'Autotrans_shift';
Br = BreachSimulinkSystem(mdl); 

% define the formula
phi1 = STL_Formula('phi1', '(alw (speed[t]<vmax)) and (alw (RPM[t]<rpm_max))');
phi1 = set_params(phi1,{'vmax', 'rpm_max'}, [160 4500]);

% define the input parameters and ranges 
input_params.names = {'throttle_u0'};
input_params.ranges = [0 100];
 
% defines the falsification problem and solve it
falsif_pb = FalsificationProblem(Br, phi1, input_params.names, input_params.ranges);

% solves using the default solver
falsif_pb.solve();

% collect the logged samples and traces 
Blog = falsif_pb.GetLog();
F = BreachSamplesPlot(Blog)
