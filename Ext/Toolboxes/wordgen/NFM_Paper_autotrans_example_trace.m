%% TA based signal generator setup 
%
init_TA_signal_gen;

S1 = S0.copy();
S1.SetParamRanges([pevts pbranching], [0 1]);
pvar = S1.GetVariables();

S1.Sim();
figure;
S1.PlotSignals();

mdl = 'Autotrans_wordgen';
Ba = BreachSimulinkSystem(mdl);
Ba.SetTime(time);
Ba.SetInputGen(S1);
Ba.SetParamRanges([pevts pbranching], [0 1]);
Ba.SampleDomain(pvar,4);
Ba.Sim();

%% Plotting
Ga = BreachSignalsPlot(Ba);

close all;
figure;
ax = subplot(4,1,1);
grid on;
Ga.AddSignals('throttle', ax, 'all');

ax = subplot(4,1,2);
grid on;
Ga.AddSignals('brake', ax, 'all');

ax = subplot(4,1,3);
grid on;
Ga.AddSignals('speed', ax, 'all');


ax = subplot(4,1,4);
grid on;
Ga.AddSignals('gear', ax, 'all');
xlabel('Time')

save2pdf('Autotrans_TA_traces.pdf');

%% 
%% Plotting
Ga = BreachSignalsPlot(Ba);

it = 2;
figure;
ax = subplot(3,1,1);
grid on;
Ga.AddSignals('throttle', ax, it );
set(gca, 'FontSize', 14, 'LineWidth',2);
set(get(gca, 'Children'), 'LineWidth',2)

ax = subplot(3,1,2);
grid on;
Ga.AddSignals('brake', ax, it);
set(gca, 'FontSize', 14, 'LineWidth',2);
set(get(gca, 'Children'), 'LineWidth',2)


ax = subplot(3,1,3);
grid on;
Ga.AddSignals('timeword', ax, it);
xlabel('Time')
set(gca, 'FontSize', 14, 'LineWidth',2);
set(get(gca, 'Children'), 'LineWidth',2)

save2pdf('Autotrans_Wordgen_Traces.pdf')