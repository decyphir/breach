%%  
%  Set some default parameters
x1_0 = 2;
x2_0 = 3; 
Mu = 1; 
Bvdp = BreachSimulinkSystem('vdp');

%% 
% Set a domain partly invalid
Bvdp.SetDomain('Mu','double', [-1 2]); % negative values for mu will lead to numerical errors
Bvdp.GridSample(4);
Bvdp.Sim();

%%
% Error message can be recovered using printStatus
Bvdp.printStatus();

%% 
% Filtering valid/invalid traces 
[Bok, Berr, Binput_err] = Bvdp.FilterTraceStatus() % Binput_err is for input generators with non satisfied contraints
                                                                                % Not relevant here.


