function [val_opt, Popt]  = SOptimProp(Sys, P, phi, opt)
%SOPTIMPROP optimizes the satisfaction of a property
% 
% Synopsis: [val_opt, Popt] = SOptimProp(Sys, P, phi, opt)
% 
% Inputs:
%  - Sys : is a system
%  - P  : is a parameter set. Parameter values in P are used for
%          initializing the optimization algorithm. It may contain many
%          parameter vectors. Trajectories don't need to be computed.
%  - phi : is a STL formula
%  - opt : is an option structure with the following fields:
%       - tspan   : the time domain computation of the trajectories. If
%                   not provided, either Sys must have a tspan field, or P
%                   must contains computed trajectories. Otherwise, an
%                   error is thrown.
%       - tau     : (Optional, default=first tspan value) time for the
%                    evaluation of phi
%       - params  : (Optional, default=P.dim) names or indexes of variable
%                   (search) parameters.
%       - lbound  : (Optional, default=P.pts-P.epsi) lower bounds for the
%                   search domain. If not provided, all parameters in
%                   params must be uncertain parameters of P.
%       - ubound  : (Optional, default=P.pts+P.epsi) upper bounds for the
%                   search domain. If not provided, all parameters in
%                   params must be uncertain parameters of P.
%       - MaxIter : Maximal number of optimization iteration.
%       - OptimType : (Optional, defaut='Max') string which indicates the
%                     type of optimization. It must be 'Max', 'Min' or
%                     'Zero'.
%       - StopWhenFound : (Optional) set to 1 to compute satisfaction for
%                         initial parameters in P0 then stops whenever a
%                         positive ('Max') or negative ('Min') solution is
%                         found.
%       - StopWhenFoundInit : same as above except that it does not
%                             necessarily compute all trajectories in P0
%       - Ninit   : (Optional, default=all parameter vectors) tries the
%                   Ninit best initial pts
%
%       - timeout : (optional), default = inf; forces convergence after 
%                   timeout s. of computation        
%    
% Outputs:
%  - val_opt : the truth value of phi for the param set Popt. It is a
%              scalar if StopWhenFound or StopWhenFoundInit is set to 1.
%              Otherwise, it is a vector of size 1 x size(P.pts,2).
% 
%  - Popt    : if StopWhenFound or StopWhenFoundInit is set to 1, and a
%              parameter vector leading to a negative (resp. positive)
%              truth value of phi is found, Popt is this parameter vector.
%              Otherwise, Popt contains the optimum found for each
%              parameter vector of P.
% 
%See also SOptimPropLog Falsify
%

if isfield(opt, 'algo')
    algo = lower(opt.algo);
else
    algo = 'nelder-mead-wait-bar';
end

switch algo
    case 'nelder-mead'
     [val_opt, Popt]  = SOptimPropNM_no_waitbar(Sys, P, phi, opt);
    case 'nelder-mead-wait-bar' % FIXME temporary     
     [val_opt, Popt]  = SOptimPropNM(Sys, P, phi, opt);
     % will implement other algo soon...    
end
