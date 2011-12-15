function index=  FindParam(Sys,param)
%  
% FINDPARAM Finds the indices of parameters given by their name for a given system
%  
% Syntax: index = FindParam(Sys, param)
%
%  
  if ~isfield(Sys,'ParamList')
    error('No parameter list ...');
  end
    
  if (~iscell(param))
     param = {param}; 
  end
  
  index = [];
  for i= 1:numel(param)
   for j = 1:numel(Sys.ParamList)
      
      test = strcmp(Sys.ParamList{j},param{i});
      if (test)
          index = [index j];
      end
    end
  end
 