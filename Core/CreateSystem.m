function Sys = CreateSystem(vars, params, p0)
% CREATESYSTEM creates a dynamicless system      


  Sys.DimX = numel(vars);
  Sys.DimU =0; 
  Sys.DimP = numel(vars)+numel(params); 
  Sys.ParamList = {vars{:} params{:}};
  Sys.x0 = zeros(1,numel(vars));  
  
  if (~exist('p0'))
      Sys.p = zeros(1, Sys.DimP);
  else
      Sys.p = p0;
  end
  
  % dynamicless system
  Sys.type = 'traces';
  
  assignin('base','Sys',Sys);
  