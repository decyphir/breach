%% Sensitivity Analysis

%% 
% <latex>\begin{frame}[fragile]{Sensitivity analysis}</latex>
% Parameter names and ranges for sensitivity computation
params = {'AF_sensor_tol','MAF_sensor_tol','fuel_inj_tol','kappa_tol','tau_ww_tol','pump_tol','Engine_Speed_u0', 'kp', 'ki'};
ranges = [ 0.99 1.01; 0.99 1.01; 0.99 1.01;0.99 1.01; 0.99 1.01;0.99 1.01; 0 1000; 0.0 0.1; 0. 0.3];
AFC_Sensi = BrAFC.copy(); AFC_Sensi.SetParamRanges(params,ranges);

% This will compute and plot sensitivities of overshoot amplitude wrt various sensor errors 
AFC_Sensi.SensiSpec(overshoot_AFC);

%% 
% <latex>\end{frame}</latex>