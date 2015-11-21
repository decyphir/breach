function b = isaSys(S)
% isaSys Returns true if the argument is a classic Breach Sys structureb

b = isstruct(S)&&all(isfield(S,{'type','DimX','DimP','ParamList','p'}));
