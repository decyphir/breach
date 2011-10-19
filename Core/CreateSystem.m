function CreateSystem(vars, params)
      
  Sys.DimX = numel(vars);
  Sys.DimU =0; 
  Sys.DimP = numel(vars)+numel(params); 
  Sys.ParamList = {vars{:} params{:}};

  Sys.x0 = zeros(1,numel(vars));
  Sys.p = zeros(1, Sys.DimP);

  Sys.type = 'traces';
  
  assignin('base','Sys',Sys);
  