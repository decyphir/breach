% In this experiment, we try to establish a template for input-output
% mapping coverage. 
% 
% The goal is to find a minimal input set that guarantees coverage of
% outputs.
%
% First we look at the reachable outputs, by bombarding with input signals,
% then we use basic filtering. 
% 

BrDemo.InitAFC

% setup input generation

B = BrAFC.copy();
B.SetParamRanges('Pedal_Angle_pulse_period',[.1 10]);
B.SetParamRanges('Pedal_Angle_pulse_amp',[0 100]);
B.SetParamRanges('Pedal_Angle_pulse_width',[0 1]);
B.SetParamRanges('Pedal_Angle_pulse_delay',[0 10]);

Engine_params = B.expand_param_name('Engine_Speed');
B.SetParamRanges(Engine_params, [900 1100]);

% Generate 
B.QuasiRandomSample(100);
B.Sim()
B.PlotSignals('AF')

%% Compute coverage
Bc = B.copy();
opts = struct;
opts = B.SetCoverageOptions(opts,'ParamsMode','off');
B.GetCoverage();

%% 
Bp = B.copy();
opc = Bp.SetCoverageOptions('ParamsMode','off', 'IncludeSignals',{'AF','MAF'});
disp_cover_opts(opc)
Bp.GetCoverage(opc)