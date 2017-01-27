function varID = MakeUniqueID(varName, varProtected)
% TODO use matlab.lang functions in replacement of genvarname for ver> R2014a

varID = genvarname(varName, varProtected);
