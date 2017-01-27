function [t, X] = sim_autotrans_accel(Sys, tspan, pts)

mdl = 'Autotrans_shift_accel';

params = pts(Sys.DimX+1:end); %  pts(1:Sys.DimX) are for signals initial conditions, irrelevant here

t1 = params(1);
throttle0 = params(2);
throttle1 = params(3);
brake0 = params(4);
brake1 = params(5);

r_inputs.time= [0 ; t1; tspan(end)];
r_inputs.signals(1).dimensions = 1;
r_inputs.signals(1).time = [ 0; t1; tspan(end)];
r_inputs.signals(1).values = [throttle0; throttle1; throttle1];
r_inputs.signals(2).dimensions = 1;
r_inputs.signals(2).time  = [0 ; t1; tspan(end)];
r_inputs.signals(2).values = [brake0; brake1; brake1];

assignin('base', 'r_inputs', r_inputs); 

r_sim_out = sim(mdl, 'LoadExternalInput', 'on',....
    'SimulationMode', 'rapid',...
    'SaveOutput', 'on',  ...
    'OutputSaveName', 'r_outputs', ...
    'ExternalInput', 'r_inputs');

out = r_sim_out.get('r_outputs');
t = out.time';
X = [ out.signals(1).values'; ...
      out.signals(2).values'; ...
      out.signals(3).values'; ...
      out.signals(4).values'; ...
      out.signals(5).values'; ...
 ];

 
