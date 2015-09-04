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

  set_param(mdl, 'OutputTimes', 'tspan',...
                 'OutputOption','SpecifiedOutputTimes'); 
  
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
  
%   
%   lg = simout.get('logsout');
%   tout = simout.get('tout')';
%   Y = simout.get('yout');    
%   X=[];
%   
%   if ~isempty(Y)
%     for i=1:numel(Y.signals)
%       xx = interp1(tout', double(Y.signals(i).values), tspan')';
%       X = [X; xx];
%     end
%   end
%   
%   for i = Sys.DimY+1:num_signals
%     sig = lg.getElement(Sys.ParamList{i});
%     xdata = interp1(sig.Values.Time',double(sig.Values.Data(:,1)),tspan, 'linear','extrap');   
%     X = [X ; xdata(1,:)]; %% FIXME: SUPPORT FOR MULTIDIMENSIONAL SIGNALS   
%   end;
%   
%   tout = tspan;
