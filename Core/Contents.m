% BREACH/CORE Routines for trajectory computation and higher-level functions
%
% Files
%   ComputeTraj       - compute trajectories for a system given initial conditions
%   ComputeTrajSensi  - compute trajectories with  sensitivities
%   CreateSystem      - creates a dynamicless system      
%   GetF              - returns the value of rhs function of the ODE dynamics
%   GetF_Traj         - returns the values of the rhs ODE function along a computed trajectory 
%   GetTrajValues     - extract values of a variable evolution in a set of trajectories 
%   SEvalProp         - Eval property for previously computed trajectories
%   SFindPropBoundary - find parameters in a set which are close to a boundary
%   SOptimProp        - optimizes the satisfaction of a property
%   sreach            - Refine a parameter set to get a good estimate of reachset with a set of trajectories 
