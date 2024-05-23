%% Init Model
init_SDT;

%% Time and requirements
time = 0:1e-10:1e-6;
STL_ReadFile('SDT_req.stl');
R =  BreachRequirement({per, notsat});

%% TA based signal generator setup 
num_evt = 70;
sg = TA_signal_gen('In1','cycle_8_12.prism',{'a','b','c','d'}, num_evt, 'pchip');

%%
S = BreachSignalGen(sg);
S.SetParam({'In1_a_val','In1_b_val','In1_c_val','In1_d_val'}, [0 1 2 1]);

%% Add input to model 
Bta = BP2.copy();
Bta.SetInputGen(S);
Bta.SetTime(time);
ts = 1e-8;
Bta.SetParam('time_scale', ts, true);
params = sg.params(1:num_evt);
ranges = repmat([0,1], numel(params),1);
Bta.SetParamRanges(params, ranges);
var = Bta.GetVariables();
Bta.SampleDomain(var, 5);
Bta.Sim();

%% Eval requirement
R.Eval(Bta);
BreachSamplesPlot(R);

%% Testing without wordgen (free dt) 
Bva = BP2.copy();
sg_va = var_cp_signal_gen('In1',num_evt+2, 'pchip');
Bva.SetInputGen(sg_va);
Bva.SetTime(time);
Bva.SetParam('time_scale', ts,true);
params_dt = Bva.expand_param_name('In1_dt.');
params_u = Bva.expand_param_name('In1_u.');
ranges_dt = repmat([0 8*ts], numel(params_dt), 1);

%
idx_u = FindParam(Bva.P, params_u);

%
p0_u = repmat([1 2 1 0]', round(71/4), 1);

%
Bva.SetParamRanges(params_dt, ranges_dt)
Bva.SetParam(params_u, p0_u)

%
Bva.SampleDomain(params_dt, 5);
Bva.Sim()


%% Plotting
Gva = BreachSignalsPlot(Bva);
Gta = BreachSignalsPlot(Bta);
close all;
figure;
ax = subplot(8,1,[1 2]);
xlabel('Random Times');
grid on;
Gva.AddSignals('In1', ax, 'all');
ax = subplot(8,1,[4 5]);
xlabel('Using Timed Words');
grid on;
Gta.AddSignals('In1', ax, 'all');
ax = subplot(8,1,[7 8]);
grid on;
Gta.AddSignals('OutSat', ax, 'all');
xlabel('\Sigma-\Delta output')
