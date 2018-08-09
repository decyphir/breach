InitBreach
%% 
% Automatic Transmission example.  
mdl = 'Autotrans_shift';

%% 
% The following creates a Breach Simulink object from model mdl with input signals parameterized
% as variable step with 2 steps.
BrSys = BreachSimulinkSystem(mdl)
BrSys.SetInputGen('VarStep2');


%% 
% Set input values (other inputs are set to 0 by default) 

% Set acceleration for 10s at 100%
SetParam(BrSys, 'throttle_u0', 100);
SetParam(BrSys, 'throttle_dt0', 10);
     
% Set brake for 20s after
SetParam(BrSys, 'brake_u0', 0);
SetParam(BrSys, 'brake_dt0', 10);
SetParam(BrSys, 'brake_u1', 325);

%% 
% Simulation and plot signals  
BrSys.SetTime(0:.01:40);
BrSys.Sim();
BrSys.PlotSignals();

%% 
% Change some values 
BrSys.SetParam({'throttle_u0','throttle_dt0'}, [10 10]);
BrSys.Sim();
BrSys.PlotSignals;

%% 
% Robust satisfaction of some STL formulas
BrSys.CheckSpec('gear[t]>0')
BrSys.CheckSpec('alw (speed[t]< 100)')
 
%% 
% Changing some parameters - will automatically simulate the model if
% necessary
BrSys.GetRobustSat('alw (speed[t]< 100)', {'throttle_u0','throttle_dt0'}, [100 30])

%% 
% PSTL formula - since throttle and dt_u0 don't change, does not recompute
% simulation
BrSys.GetRobustSat('ev_[0, tau] (RPM[t]> pi)', {'throttle_u0','throttle_dt0', 'tau','pi'}, [100 30 20 4000])

%% 
% Getting a direct mapping from param values to robust satisfaction
phi= STL_Formula('phi','ev_[0, tau] (RPM[t]> pi)'); %slightly better to define phi before
robfn = BrSys.GetRobustSatFn(phi,{'throttle_u0','throttle_dt0', 'tau','pi'}); 

%% 
x= [100 30 20 4000]
robfn(x) 
% same result as before

%% 
x= [100 30 20 3000]
robfn(x)
