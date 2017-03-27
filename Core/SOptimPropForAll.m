function [val_opt, Pout, timed_out]  = SOptimProp(Sys, P, phi, opt)
%
% SOPTIMPROPFORALL maximizes the satisfaction of a property for a set of parameters
%
% Synopsis: [val_opt, Popt]  = SOptimPropForAll(Sys, P0, phi, opt)
%
% Input:
%    - Sys is a system
%    - P0 is a parameter set for Sys. Parameter values in P0 are used for
%                  initializing the optimization algorithm
%    - phi is a STL property
%    - opt is an option structure with the following fields :
%
%        - tspan   : the time domain computation of the trajectories. If
%                    not provided, either Sys must have a tspan field, or P
%                    must contains computed trajectories. Otherwise, an
%                    error is thrown.
%        - tau     : time for the evaluation of phi (default = first tspan
%                    value)
%        - params  : variable (search) parameters. If not provided, params
%                    is based on P.dim.
%        - lbound  : lower bounds for the search domain. If not provided,
%                    all parameters in params must be uncertain parameters
%                    of P, and lbound is defined as P.pts-P.epsi.
%        - ubound  : upper bounds for the search domain. If not provided,
%                    all parameters in params must be uncertain parameters
%                    of P, ans ubound is defined as P.pts+P.epsi.
%        - MaxIter : (mandatory) max number of optimization iteration. If
%                    not provided, an error is thrown.
%        - StopWhenFound : set to 1 to compute satisfaction for initial
%                          parameters in P0 then stops whenever  a positive
%                          ('Max') or negative ('Min') solution is found
%        - Ninit   : tries the Ninit best initial pts
%
%
% Output:
%    - val_opt : the truth value of phi for the param set Sopt. It is a
%                scalar if StopWhenFound or StopWhenFoundInit is set to 1.
%                Otherwise, it is a vector of size 1 x size(P.pts,2).
%
%    - Popt    : if StopWhenFound or StopWhenFoundInit is set to 1, and a
%                set of parameter values leading to a negative (resp.
%                positive) truth value of phi is found, Popt is this
%                parameter set.
%                Otherwise, it contains the optimum found for each set of
%                parameter values in P.
%
% See also SOptimPropLog Falsify
%

tic;

%% process options

global Ptmp; % temporary param set used to get non-variables parameter values in optim func
global found; % non empty if we found a positive or negative truth value of prop
global StopWhenFound; % cf doc
global fopt; % best truth value of prop found for the current initial set of value in P
global xopt; % parameter vector leading to the optimal satisfaction
global dim;
global timeout;
global Popt;

timeout=0;

if isfield(opt, 'tspan')
    tspan = opt.tspan;
elseif isfield(Sys, 'tspan')
    tspan = Sys.tspan;
elseif isfield(P, 'traj')
    tspan = P.traj{1}.time;
else
    error('SOptimeProp:noTspan','The field opt.tspan is not provided.');
end

if isfield(opt,'tau')
    tau = opt.tau;
    if tspan(1) > tau
        tspan = [tau tspan];
    end
else
    tau = tspan(1);
end

if isfield(opt,'params')
    dim = FindParam(P,opt.params);
    dim = dim(dim<size(P.pts,1)); % keep only existing parameters (either system or constraint parameter)
    P0 = CreateParamSet(Sys, dim);
    P0 = QuasiRefine(P0, size(P.pts,2));
    P0.pts = P.pts;
    P  = P0;
else
    dim = P.dim;
end

if isfield(opt,'ubound')
    ubound = opt.ubound;    
else
    ubound = [];
end

if isfield(opt,'lbound')
    lbound = opt.lbound;    
else
    lbound=[];
end

if isfield(opt,'MaxIter')
    MaxIter = opt.MaxIter;
else
    error('SOptimProp:noMaxIter','The field opt.MaxIter is not provided.');
end

if isfield(opt,'StopWhenFound')
    StopWhenFound = opt.StopWhenFound;
else
    StopWhenFound = 0;
end

if isfield(opt,'timeout')
    max_time=opt.timeout;
else
    max_time=inf;
end

phi = STL_OptimizePredicates(Sys,phi); % optimization of the predicates
fun = @(x) fun_rob(x, max_time, Sys, phi, tspan, tau);

%% Initial values

Ps = PdecomposeByDim(P, dim);
nb_init = numel(Ps);
found = [];
fopt = -inf;
xopt= [];

rfprintf_reset();
val = zeros(1,nb_init);
for ii = 1:nb_init
    Ptmp = Ps(ii) ;
    x = GetParam(Ptmp, dim);
    x = x(:,1);
    val(ii) = fun(x);
    
    status = ['Init ' num2str(ii) '/' num2str(nb_init) ' Robustness value: ' num2str(-val(ii))];
    rfprintf(status);
    
    if ((val(ii)<0)&&(StopWhenFound))
        val_opt = -val(ii);
        Pout = Ptmp;
        fprintf('\nSatisfying parameter vector found.\n');
        timed_out = 0;
        return;
    end
end

[~, iv] = sort(val);

%% Main Loop

options = optimset('MaxIter', MaxIter);
kk=0;
Popt = Ps(iv(1));

for ii = iv(1:nb_init)
    kk = kk+1;
    
    fprintf('\nOptimize from init point %d/%d Initial value: %g\n', kk, numel(iv), -val(ii));
    rfprintf_reset();
    
    Ptmp = Ps(ii);
    x0 = GetParam(Ptmp, dim);
    x0 = x0(:,1);
    
    [~, fopt] = optimize(fun,x0,lbound,ubound,[],[],[],[],[],[],options,'NelderMead');
    fprintf('\n');
    
    if((StopWhenFound)&&(~isempty(found)))
        break;
    end
end

Pout = Popt; 
timed_out = timeout;
if isfield(Sys, 'init_fun') % init_fun can modify non-uncertain parameters
    Popt = Sys.init_fun(Popt);
end
if isfield(Popt, 'init_fun')
    Popt = Popt.init_fun(Popt);
end
val_opt=-fopt;

end

function val = fun_rob(x, max_time, Sys, phi, tspan, tau)
%% function fun_rob
global Ptmp fopt Popt found StopWhenFound timeout xopt dim

if(StopWhenFound&&~isempty(found)) % positive value found, do not need to continue
    val = -found; % optimize tries to minimize the objective function, so
    return ;          % we provide it the opposite of the truth value
end

ct= toc;
if (ct>max_time)
    timeout=1;
    val = -fopt;     % forces convergence of the optimizer in case timeout occured;
    return;
end

Ptmp= SetParam(Ptmp,dim, x);
Ptmp = ComputeTraj(Sys, Ptmp, tspan);

val = min(STL_Eval(Sys, phi, Ptmp, Ptmp.traj, tau));

if(val>0)
    found = val;
    fprintf('Found !\n');
end

if(val>fopt)
    fopt = val;
    Popt = Ptmp;
    xopt = GetParam(Ptmp,dim); % as ComputeTraj launch init_fun, Ptmp.pts can be different than x
    xopt = xopt(:,1);
end

status = ['Robustness value: ' num2str(val) ' Current optimal: ' num2str(fopt) ' Computation Time: ' num2str(ct)];
fprintf(status);
fprintf('\n');
val = -val; % optimize tries to minimize the objective function, so we
% provide it -val instead of val
end
