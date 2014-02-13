%CVM is a function to call CVODES functions.
% 
% Synopsis: P = cvm(mode1, P, tspan, u)
%       OR:     cvm(mode2, Ns, meth, yS0, options)
% 
% Inputs:
%  - mode1 : mode is an integer describing the called function. Possible
%            values for mode are 61 (compute trajectories), 93 (compute
%            trajectories and sensitivities)
%  - mode2 : 2 (initialisation), 
%  - P     : 
%  - tspan : 
%  - u     : 
% 
% Output:
%  - P : parameter set augmented with corresponding fields
% 
% Credits:
% -- Radu Serban @ LLNL -- April 2005
% 
%See also ComputeTraj ComputeTrajSensi
