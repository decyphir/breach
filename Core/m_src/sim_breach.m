function [tout, X] = sim_breach(Sys, tspan, pts)
% 
% Generic wrapper function that runs a Simulink model and collect signal
% data in Breach format (called by ComputeTraj) 
%  

  mdl = Sys.mdl;
  load_system(mdl);
  num_signals = Sys.DimX;
  
  params = Sys.ParamList;   
  for i = 1:numel(params)-num_signals
    assignin('base',params{i+num_signals},pts(i+num_signals));
  end
  
  assignin('base','tspan',tspan);

 if numel(tspan)>2 
  set_param(mdl, 'OutputTimes', 'tspan',...
                 'OutputOption','SpecifiedOutputTimes'); 
 else
  set_param(mdl, 'OutputTimes', 'tspan',...
                 'OutputOption','RefineOutput');      
 end 

  try                
    simout= sim(mdl);
  catch
    s= lasterror;
    warning(['An error was returned from Simulink:' s.message '\n Returning a null trajectory']);
    tout = tspan;
    X = zeros(Sys.DimX, numel(tspan)); 
    return;
  end
  
  [tout, X] = simout2X(simout);