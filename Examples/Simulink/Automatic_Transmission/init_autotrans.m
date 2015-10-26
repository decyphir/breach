% Init Model with two inputs 

mdl = 'Autotrans_shift';
Sys = CreateSimulinkSystem(mdl, 'logged', {}, [], 'VarStep2');
Sys.tspan = 0:.1:40; %default time span

% Set default input values (other than 0) 

% Accelerate for 20 s at 100%
Sys = SetParam(Sys, 'throttle_dt0', 20);
Sys = SetParam(Sys, 'throttle_u0', 100);
Sys = SetParam(Sys, 'brake_dt0', 20);
Sys = SetParam(Sys, 'brake_u0', 0);
     
% Brake for 20s after
Sys = SetParam(Sys, 'throttle_dt1', 20);
Sys = SetParam(Sys, 'throttle_u1',   0);
Sys = SetParam(Sys, 'brake_dt0', 20);
Sys = SetParam(Sys, 'brake_u0', 325);
