% Init Model with two inputs 

mdl = 'Autotrans_shift';
Sys = CreateSimulinkSystem(mdl, 'logged', {}, [], 'VarStep2');
Sys.tspan = 0:.1:50; %default time spa

% Set default input values (other than 0) 

% Accelerate for 20 s at 100%
Sys = SetParam(Sys, 'dt_u0', 20);
Sys = SetParam(Sys, 'throttle_u0', 100);
     
% Brake for 20s after
Sys = SetParam(Sys, 'dt_u1', 20);
Sys = SetParam(Sys, 'brake_u1', 325);

