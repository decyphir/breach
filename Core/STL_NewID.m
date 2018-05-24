function id = STL_NewID(id)

InitBreach
global BreachGlobOpt;

if isfield(BreachGlobOpt, 'STLDB')
    id  = matlab.lang.makeUniqueStrings(id, BreachGlobOpt.STLDB.keys);
end
    
    