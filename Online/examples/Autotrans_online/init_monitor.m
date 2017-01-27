mdl = 'Autotrans_online';

% Parameters
max_rob = 10000;              % initial robustness is between -max_rob, +max_rob
rob_up_lim = -max_rob;  % stops when rob_up is below 
rob_low_lim = max_rob; 
rob_diff_lim = -inf;

phi_autotrans = 'ev_[0,50] gear[t]>200';

AT = BreachSimulinkSystem(mdl);

AT.SetParam('Throttle_u0', 100);

AT.Sim();
AT.PlotSignals()





