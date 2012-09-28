function [tout X] = sim_breach(Sys, tspan, pts)
%
%
%

  mdl = Sys.mdl;
  num_signals = Sys.DimX;  
  params = Sys.ParamList;
     
  for i = 1:numel(params)-num_signals
    assignin('base',params{i+num_signals},pts(i+num_signals));
  end
  
  assignin('base','tspan',tspan);

  set_param(mdl, 'OutputTimes', 'tspan',...
                     'OutputOption', 'SpecifiedOutputTimes'); 
  
  simout= sim(mdl);              

  lg = simout.get('logsout');
  
  if ~isempty(lg)
    lg.unpack('all');
  end
  
  tout = simout.get('tout')';
  X = simout.get('yout')';    
  if ~isempty(X)
    X = [X.signals.values]';
  end
  
  for i = Sys.DimY+1:num_signals
    sig = Sys.ParamList{i};
    xdata = eval([sig '.Data']);
    xtime = eval([sig '.Time']);
    X = [X; interp1(xtime',double(xdata(:,1)'),tout)];
  end;  
