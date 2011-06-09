function index=  FindParam(Sys,param)
%  
% Finds the indices of parameters given by their name for a given system
%  
% Syntax: index = FindParam(Sys, param)
%
%  
  if ~isfield(Sys,'ParamList')
    error('No parameter list ...');
  end
    
  if (iscell(param))
     param = char(param); 
  end
  for j = 1:numel(Sys.ParamList)
      
      test = strcmp(Sys.ParamList{j},param);
      if (test)
          index = j;
          return;
      end
  end
  
  index =0;
  