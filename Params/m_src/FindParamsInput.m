function params_u = FindParamsInput(Sys)

InputNames = Sys.ParamList(Sys.DimX-Sys.DimU+1:Sys.DimX);
U = Sys.init_u(InputNames, [], []);
params_u= U.params;