
% Model
mdl = 'Autotrans_shift';
Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'UniStep1');
Sys.tspan = 0:.01:50;

%formulas = QMITL_ReadFile('spec.stl');
QMITL_Formula('phi1', '(alw (speed[t]<vmax)) and (alw (RPM[t]<rpm_max))');

params_prop.names = {'vmax', 'rpm_max'};
params_prop.values = [160 5000];

falsif_opt.params = {'throttle_u0'};
falsif_opt.ranges = [0 100];
  
falsif_opt.nb_init = 10;
falsif_opt.nb_iter = 100;

Pf = Falsify(Sys, phi1, falsif_opt, params_prop);