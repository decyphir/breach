%% testing corners with events
%
num_evt = 5;
init_TA_signal_gen;

S1 = S0.copy();
S1.SetParamRanges([pevts pbranching], [0.01 0.99]);
pvar = S1.GetVariables();
S1.CornerSample(2);

S1.Sim();
figure;
S1.PlotSignals();
