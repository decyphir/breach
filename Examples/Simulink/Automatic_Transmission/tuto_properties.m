%% In this script, we demonstrate the monitoring of a property

%% Init system and nominal trajectory (in P0)

init_autotrans;
P0 = CreateParamSet(Sys); % nominal parameter set for Sys
P0 =  ComputeTraj(Sys, P0, Sys.tspan)

%%  define property phi1

QMITL_Formula('phi1', '(alw (speed[t]<vmax)) and (alw (RPM[t]<rpm_max))');

%% Eval and plot property phi1

% phi1 has two parameters vmax and rpm_max, so we need to define them for
% P0
params_phi.names = {'vmax', 'rpm_max'};
params_phi.values =  [140 4000];

P0 = SetParam(P0, params_phi.names, params_phi.values);

% Eval satisfaction function for phi1, P0 and P0.traj at time tau=0
tau=0;
val = QMITL_Eval(Sys, phi1, P0, P0.traj, tau);

% computes and plots satisfaction function for all tau>0
figure;
SplotSat(Sys, P0, phi1 );
pause  

% plot satisfaction function with sub formula 
figure;
depth = 3;
SplotSat(Sys, P0, phi1, depth);
pause

% plot satisfaction function as a function of the system input
Pu = CreateParamSet(Sys, {'dt_u0'}, [0 10], 100); % 100 samples with dt_u0 in range [0 10]
Pu = ComputeTraj(Sys, Pu, Sys.tspan);   % Pu.traj now contains 100 trajectories 

Pu = SetParam(Pu, params_phi.names, params_phi.values); % set property parameters vmax and rpm_max
[Pu, valu] = SEvalProp(Sys, Pu, phi1);  % SEvalProp calls QMITL_Eval for all trajectories in Pu
                                        % valu(i) contains the satisfaction
                                        % value of traj(i)      
figure;
SplotProp(Pu, phi1); 



