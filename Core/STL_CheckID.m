function [b, phi] = STL_CheckID(id)
% 0 doesn't exist, 1 predicate, 2 formula, 3 explicit predicate (generated id)
b=0; 
phi = [];
InitBreach
global BreachGlobOpt;

if isfield(BreachGlobOpt, 'STLDB')
    if BreachGlobOpt.STLDB.isKey(id)
       if strcmp(get_type(BreachGlobOpt.STLDB(id)),'predicate')
           if ~isempty(regexp(id,'.+__$', 'once'))
               b=3;
           else
               b=1;
           end
       else
           b=2;
       end
        phi = BreachGlobOpt.STLDB(id);
    end
end
    
    