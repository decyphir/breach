function b = isaP(S)
% isaP Returns true if the argument is a classic Breach param set structure

b = isstruct(S)&&all(isfield(S,{'DimX','DimP','ParamList','pts'}));
