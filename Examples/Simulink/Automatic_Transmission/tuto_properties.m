%% In this script, we demonstrate the monitoring of a property

%% Init system (creates object Br) 
init_autotrans;

%% Computes and plots nominal trajectory
Br.Sim()
Br.PlotSignals()

%%  define property phi1
phi1 = STL_Formula('phi1', '(alw (speed[t]<vmax)) and (alw (RPM[t]<rpm_max))');


%% Eval and plot property phi1

% phi1 has two parameters vmax and rpm_max, so we need to define them 
phi1 = set_params(phi1,{'vmax', 'rpm_max'}, [140 4000]);

% Eval satisfaction function for phi1, P0 and P0.traj at time tau=0
val = Br.CheckSpec(phi1);

% computes and plots satisfaction function for all tau>0
Br.PlotRobustSat(phi1)
 
% plot satisfaction function as a function of system input
Br.PlotRobustMap(phi1, ...
                 {'throttle_dt0'},... % input parameter to vary
                 [0 10], ...          % input range 
                 50);                 % number of input values tested 

% plot 2d maps
Br.PlotRobustMap(phi1, ...
                 {'throttle_u0', 'throttle_dt0'},... % input parameter to vary
                 [0 100; 0 10], ...                  % input range 
                 [10 10]);                           % 15x15 grid 
            