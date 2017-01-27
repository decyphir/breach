function st = SysName(Sys)
%
% SysName   get system name based on its folder
%
% Synopsis:  st = SysName(Sys)
%  
%    

if isfield(Sys,'mdl')
  st = Sys.mdl;
else
  dr = Sys.Dir;
  indst = strfind(dr, filesep);
  st = dr(indst(end)+1:end);    
end
  
