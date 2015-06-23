

mdl = 'helicopter';

K = 10;
i = 0; 
a = 1;


Sys = CreateSimulinkSystem(mdl)
Sys.tspan = 0:.01:4;
Sys = SetParam(Sys, 'psi_u0', 10);