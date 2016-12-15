% Run Init_AFC_online first
phi_AFC = 'alw_[10, 30] ( ( (abs(AF[t]-AFref[t]) > 0.1)) => (ev_[0, 5] (abs(AF[t]-AFref[t]) < 0.1)))';

BrAFC.ResetSimulations();
BrAFC.SetParam({'max_rob','Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [.5, 12, 50]);
BrAFC.Sim(0:.1:40) 

Trace = BrAFC.GetTraces();
idx = FindParam(BrAFC.Sys, {'rob_low', 'rob_up'});
t = Trace.time;

%% 
close 
figure;
subplot(3,1,1)

plot(t,Trace.X(1,:)',t(2:end),Trace.X(2,2:end)', 'LineWidth', 4);
set(gca, 'LineWidth', 2, 'FontSize',18)

legend({'AF','AFref'});
grid on;

g = title(phi_AFC);
set(g,'Interpreter','None')
subplot(3,1,[2 3]);
hold on;
stairs(t, Trace.X(idx(2),:)', 'g', 'LineWidth', 2);
stairs(t, Trace.X(idx(1),:)', 'r', 'LineWidth', 2);
set(gca, 'LineWidth', 2, 'FontSize',14)
legend({'Upper robustness','Lower robustness'});
plot(t,0*t,'k')
grid on;

%%
%save2pdf('RobustOnlinePlot.pdf')    