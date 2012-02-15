function st = SysName(Sys)
%
% SysName   get system name based on its folder
%
% Synopsis:  st = SysName(Sys)
%  
%    
  
  dr = Sys.Dir; 
  indst = strfind(dr, filesep);
  st = dr(indst(end)+1:end);    

  