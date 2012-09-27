function index =  FindParam(Sys,param)
%  
% FINDPARAM Finds the indices of parameters given by their name for a given system or param set
%  
% Syntax: index = FindParam(Sys, param)
%
% param can be a string or a cell of string when looking for more than one
% parameter. Returns index greater than Sys.DimP for parameter(s) not found.
%  
% Example: 
%  
%  CreateSystem;
%  idx = FindParam(Sys, 'a');
%  idxs = FindParam(Sys, {'a','F','G'}); 
%
  
  if ~isfield(Sys,'ParamList')
    error('No parameter list ...');
  end
    
  if (~iscell(param))
     param = {param}; 
  end
  
  index = ones(1, numel(param));

  newp = Sys.DimP;
  
  if isfield(Sys, 'pts')
      newp= size(Sys.pts,1);
  end
  for i= 1:numel(param)
    found = 0;
    for j = 1:numel(Sys.ParamList)
     
     test = strcmp(Sys.ParamList{j},param{i});
     if (test)
       found = 1; 
       index(i) = j;
     end
   end
   if (~found)
     newp = newp+1;
     index(i) = newp;
   end     
  end
 
  