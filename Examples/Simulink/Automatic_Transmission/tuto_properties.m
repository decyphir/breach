%% Init system and nominal trajectory (in P0)

tuto_simulink_system;
close;

%%  define property phi1

QMITL_Formula('phi1', '(alw (speed[t]<vmax)) and (alw (RPM[t]<rpm_max))');

%% Eval and plot property phi1

% phi1 has two parameters vmax and rpm_max, so we need to define them for
% P0

P0 = SetParam(P0, {'vmax', 'rpm_max'}, [140 4000]);


% Eval satisfaction function for phi1, P0 and P0.traj at time tau=0

val = QMITL_Eval(Sys, phi1, P0, P0.traj);

% plot satisfaction function for all tau>0

figure;
SplotSat(Sys, P0, phi1 );
  
% plot satisfaction function with sub formula 

figure;
depth = 3;
SplotSat(Sys, P0, phi1, depth);

%  BreachGlobOpt.RobustSemantics = 1