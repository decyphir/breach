% New copy of the interface object and specification to falsify
AFC_Optimize = BrAFC.copy();

% We specify the parameters of the specification to falsify in a structure:
params_prop.names = {'ti', 'tf', 'tau', 'c'};
params_prop.values = [  10, 40, 3 ,  .01*14.7];

AFC_Optimize.SetParam(params_prop.names,params_prop.values);

problem = BreachProblem();
problem.options = optimset('Display', 'iter', 'TolFun', 1e-4);

% Next we specify the parameters and ranges for the system.   
params = {'Pedal_Angle_base_value','Pedal_Angle_pulse_period'};
AFC_Optimize.SetParamRanges(params, [0 61.1; 10 15]);

AFC_Optimize.MaxSatSpec(overshoot_AFC, params, problem);

