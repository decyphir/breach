%% TA based signal generator setup 
%
init_TA_signal_gen;

S1 = S0.copy();
S1.SetParamRanges([pevts pbranching], [0 1]);
pvar = S1.GetVariables();
%S1.SampleDomain(pvar,100);

S1.Sim();
figure;
S1.PlotSignals();

%%  Checking time duration

dt_evt = S1.GetParam(pevts);

%%
T_evt = sum(dt_evt, 1);
mean(T_evt)