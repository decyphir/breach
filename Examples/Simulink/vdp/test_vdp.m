%% 
x1_0 = 2;
x2_0 = 3; 
Mu = 1; 
Bvdp = BreachSimulinkSystem('vdp');

%% 
Bvdp.SetDomain('Mu','double', [-1 2]);
Bvdp.GridSample(4);
Bvdp.Sim();
[Bok, Berr, Binput_err] = Bvdp.FilterTraceStatus();

