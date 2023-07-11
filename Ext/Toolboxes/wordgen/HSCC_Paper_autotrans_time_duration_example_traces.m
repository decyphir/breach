%% TA based signal generator setup 
%
init_TA_signal_gen;

S1 = S0.copy();
S1.SetParamRanges([pevts pbranching], [0 1]);
pvar = S1.GetVariables();

S1.Sim();
figure;
S1.PlotSignals();

mdl = 'Autotrans_shift';
Ba = BreachSimulinkSystem(mdl);
Ba.SetTime(time);
Ba.SetInputGen(S1);
Ba.SetParamRanges([pevts pbranching], [0 1]);
Ba.SampleDomain(pvar,10);
Ba.Sys.Verbose=0;
%%
Ba.Sim();

%% Plotting
Ga = BreachSignalsPlot(Ba);

%%
Ga.AddSignals('throttle', 1, 'all');

%%

Ga.AddAxes
Ga.AddSignals('brake', 2, 'all');

Ga.AddAxes
Ga.AddSignals('speed', 3, 'all');

Ga.AddAxes
Ga.AddSignals('gear', 4, 'all');
xlabel('Time')

%save2pdf('Autotrans_TA_traces.pdf');

% alw (inSafeRegion[t] == 1)