function b = isaSys(S)
% isaSys Returns true if the argument is a classic Breach Sys structureb

b = isstruct(S)&&all(isfield(S,{'DimX','DimP','ParamList','p','type','x0'}));
