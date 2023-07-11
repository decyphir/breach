InitBreach

%% TA based signal generator setup 
num_evt = 10;
%sg2 = TA_signal_gen2({'throttle', 'brake'},'driving_NFM_v2.prism',{'a','b','c','d','e','f','g','h'}, num_evt);
sg2 = TA_signal_gen2({'throttle', 'brake'},'driving_HSCC21_v1.prism',{'a','b','c','d','e','f','g','h'}, num_evt);
sg2.verbose = true;
S0 = BreachSignalGen(sg2);
ts = 1;
time = 0:.005:100;
sg2.min_dt = .01;

S0.SetParam('time_scale',ts);
S0.SetTime(time);

%%  Assign values
p_throttle_init = S0.expand_param_name('.*throttle_init_val');
p_brake_init = S0.expand_param_name('.*brake_init_val');

p_throttle_a = S0.expand_param_name('.*throttle_a_val');
p_brake_a = S0.expand_param_name('.*brake_a_val');

p_throttle_b = S0.expand_param_name('.*throttle_b_val');
p_brake_b = S0.expand_param_name('.*brake_b_val');

p_throttle_c = S0.expand_param_name('.*throttle_c_val');
p_brake_c = S0.expand_param_name('.*brake_c_val');

p_throttle_d = S0.expand_param_name('.*throttle_d_val');
p_brake_d = S0.expand_param_name('.*brake_d_val');

p_throttle_e = S0.expand_param_name('.*throttle_e_val');
p_brake_e = S0.expand_param_name('.*brake_e_val');

p_throttle_f = S0.expand_param_name('.*throttle_f_val');
p_brake_f = S0.expand_param_name('.*brake_f_val');

p_throttle_g = S0.expand_param_name('.*throttle_g_val');
p_brake_g = S0.expand_param_name('.*brake_g_val');

p_throttle_h = S0.expand_param_name('.*throttle_h_val');
p_brake_h = S0.expand_param_name('.*brake_h_val');


%%
acc_range = [0 100];
brake_range = [100 325];

%%
S0.SetParamRanges(p_throttle_init, acc_range);
S0.SetParam(p_brake_init, 0);

%% Ranges for label a 
S0.SetParam(p_throttle_a, 0);
S0.SetParam(p_brake_a, 0);

%% Ranges for label b 
S0.SetParam(p_throttle_b, 0);
S0.SetParamRanges(p_brake_b, brake_range);

%%
S0.SetParamRanges(p_throttle_c, acc_range);
S0.SetParam(p_brake_c, 0);

S0.SetParam(p_throttle_d, 0);
S0.SetParamRanges(p_brake_d, brake_range);

S0.SetParam(p_throttle_e, 0);
S0.SetParamRanges(p_brake_e, brake_range);

S0.SetParam(p_throttle_f, 0);
S0.SetParam(p_brake_f, 0);

S0.SetParamRanges(p_throttle_g, acc_range);
S0.SetParam(p_brake_g, 0);

S0.SetParamRanges(p_throttle_h, acc_range);
S0.SetParam(p_brake_h, 0);

pvar = S0.GetVariables();
pevts = S0.expand_param_name('e.*_dt');
pbranching = S0.expand_param_name('e.*_branching');

