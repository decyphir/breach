function phi = set_id(phi,id)
% SET_ID assign a new id to a formula

global BreachGlobOpt
phi.id=id;
BreachGlobOpt.STLDB(phi.id) = phi;