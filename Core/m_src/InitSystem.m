function InitSystem(Sys)
%INTISYSTEM initialises the system with t0 = 0, x0 = Sys.x0, ... This
% function calls CVodeMalloc.
% 
% Synopsis: InitSystem(Sys)
%

t0 = 0.;
if isfield(Sys,'fdata')
    data = Sys.fdata;
end

data.p = Sys.p;
options = Sys.CVodesOptions;
x0 = Sys.x0;

CVodeMalloc([],t0,x0,options,data);

end
