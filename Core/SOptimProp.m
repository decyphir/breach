function [val_opt, Popt]  = SOptimProp(Sys, P, phi, opt)
%
% SOPTIMPROP optimizes the satisfaction of a property
%
% Synopsis: [val_opt, Popt]  = SOptimProp(Sys, P0, phi, opt)
%
% Input:
%    - P0 is a parameter set for Sys. Parameter values in P0 are used for
%                  initializing the optimization algorithm
%    - phi is a QMITL property
%    - opt is an option structure with the following fields :
%
%        - tspan  : the time domain computation of the trajectories
%        - tau    : time for the evaluation of phi (default first tspan
%                   value)
%        - params : (mandatory) variable (search) parameters
%        - lbound : (mandatory) lower bounds for the search domain
%        - ubound : (mandatory) upper bounds for the search domain
%        - MaxIter : (mandatory) max number of optimization iteration
%        - OptimType : 'Max' (default), 'Min' or 'Zero'
%        - StopWhenFound : set to 1 to compute satisfaction for initial
%                          parameters in P0 then stops whenever  a positive
%                          ('Max') or negative ('Min') solution is found
%        - StopWhenFoundInit : same as above except that it does not
%                              necessarily compute all trajectories in P0
%        - Ninit : tries the Ninit best initial pts
%
%
% Output:
%    - val_opt : the truth value of phi for the param set Sopt. It is a
%                scalar if StopWhenFound or StopWhenFoundInit it sets to 1.
%                Otherwize, it is a vector of size 1 x size(P.pts,2).
%
%    - Popt    : if StopWhenFound or StopWhenFoundInit is set to 1, and a
%                set of parameter values leading to a negative (resp.
%                positive) truth value of phi is found, Popt is this
%                parameter set.
%                Otherwize, it contains the optimum found for each set of
%                parameter values in P.
%
% See also SOptimPropLog
%

%% process options

global Ptmp; % temporary param set used to get non variables parameter values in optim func
global found; % non empty if we found a positive or negative truth value of prop
global StopWhenFound; % cf doc
global fopt; % best truth value of prop found for the current initial set of value in P
global traj_opt; % trajectory leading to fopt truth value of prop
global xopt; % parameter set leading to the trajectory traj_opt

found = [];
if isfield(opt, 'tspan')
    tspan = opt.tspan;
elseif isfield(Sys, 'tspan')
    tspan = Sys.tspan;
elseif isfield(P, 'traj')
    tspan = P.traj(1).time;
else
    tspan = 0:.2:10;
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
else
    dim = P.dim;
end

if isfield(opt,'OptimType')
    OptimType = lower(opt.OptimType); % to avoid case mistake, we convert to lower case
else
    OptimType = 'max';
end

if isfield(opt,'MaxIter')
    MaxIter = opt.MaxIter;
end

if isfield(opt,'StopWhenFound')
    StopWhenFound = opt.StopWhenFound;
else
    StopWhenFound = 0;
end

if isfield(opt,'StopWhenFoundInit')
    StopWhenFoundInit = opt.StopWhenFoundInit;
    StopWhenFound = StopWhenFound | StopWhenFoundInit;
else
    StopWhenFoundInit = 0;
end

if isfield(opt,'Ninit')
    Ninit = opt.Ninit;
else
    Ninit = size(P.pts, 2);
end


%% Initial values
options = optimset('MaxIter', MaxIter);

if (StopWhenFoundInit)
    nb_errors = 0; % number of ComputeTraj errors
    rfprintf_reset();
    val = zeros(1,size(P.pts,2));
    for i = 1:size(P.pts, 2)
        Ptmp = Sselect(P,i);
        try
            Ptmp = ComputeTraj(Sys, Ptmp, tspan);
            val(i) = QMITL_Eval(Sys, phi, Ptmp, Ptmp.traj, tau);
        catch % in case an error occurs during computation of ComputeTraj
            warning('SOptimProp:ComputeTraj','Error during computation of an initial trajectory, keep going.');
            if strcmp(OptimType,'max')
                val(i) = -inf;
            else
                val(i) = inf;
            end
            nb_errors = nb_errors + 1;
        end
        status = ['Init ' num2str(i) '/' num2str(size(P.pts, 2)) ' Robustness value: ' num2str(val(i))];
        rfprintf(status);
        
        switch OptimType
            case 'max'
                if val(i)>0
                    val_opt = val(i);
                    Popt = Ptmp;
                    fprintf('\n'); % to have a pretty display
                    return;
                end
            case 'min'
                if val(i)<0
                    val_opt = val(i);
                    Popt = Ptmp;
                    fprintf('\n'); % to have a pretty display
                    return;
                end
        end
        
        if i==1
            Popt = Ptmp;
        else
            Popt = SConcat(Popt, Ptmp);
        end
    end
    Ninit = min(Ninit,size(P.pts,2)-nb_errors); % we dont consider parameter sets generating an error of ComputeTraj
else
    try
        Popt = ComputeTraj(Sys, P, tspan);
    catch err
        fprintf([err.message,'\n']);
        error('SOptimProp:ComputeTraj','Error during computation of initial trajectories. Try with opt.StopWhenFoundInit=1.')
    end
    [Popt, val] = SEvalProp(Sys, Popt, phi, tau);
    Ptmp = Sselect(Popt,1);
end


switch OptimType
    case 'max'
        [~, iv] = sort(-val);
        if val(iv(1))>0 % if the highest value is positive
            found = val(iv(1));
        end
        fun = @(x) fun_max(x, Sys, phi, tspan, tau);
        
    case 'min'
        [~, iv] = sort(val);
        if val(iv(1))<0 % if the lowest value is negative
            found = val(iv(1));
        end
        fun = @(x) fun_min(x, Sys, phi, tspan, tau);
        
    case 'zero'
        [~, iv] = sort(abs(val));
        fun = @(x) fun_zero(x, Sys, phi, tspan, tau);
end

if ((StopWhenFound)&&(~isempty(found))) || (MaxIter==0)
    Popt = Sselect(Popt,iv(1));
    val_opt = found;
    return ;
end

%% Main Loop

Ninit = min(Ninit,numel(iv));
val_opt = zeros(1,Ninit);
k=0;
for i = iv(1:Ninit)
    k = k+1;
    if isfield(opt, 'lbound')
        lbound = opt.lbound;
    else
        lbound = P.pts(dim,i)-P.epsi(:,i);  % ERROR : THE ORDER OF params AND P.dim MAY DIFFER !!
    end
    
    if isfield(opt, 'ubound')
        ubound = opt.ubound;
    else
        ubound = P.pts(dim,i)+P.epsi(:,i);  % ERROR : THE ORDER OF params AND P.dim MAY DIFFER !!
    end
    
    fprintf('\nOptimize from init point %d/%d Initial value: %g\n', k, numel(iv), val(i));
    rfprintf_reset();
    x0 = Popt.pts(dim,i);
    fopt = val(i); % we initialize with the only truth value computed for this set of values
    traj_opt = Popt.traj(Popt.traj_ref(i));           % <--- !!! NOT SURE OF THAT (but I guess it is correct)
    xopt = Popt.pts(dim,i);
    [~, val_opt(k)] = optimize(fun,x0,lbound,ubound,[],[],[],[],[],[],options,'NelderMead');
    fprintf('\n');
    Popt.pts(dim,i) = xopt;
    Popt.traj(Popt.traj_ref(i)) = traj_opt;
    Popt.Xf(:,i) = traj_opt.X(:,end);
    
    if (StopWhenFound)&&(~isempty(found))
        Popt = Sselect(Popt,i);
        if isfield(Sys, 'init_fun') % init_fun can modify non-uncertain parameters
            Popt = Sys.init_fun(Popt);
        end
        if isfield(Popt, 'init_fun')
            Popt = Popt.init_fun(Popt);
        end
        val_opt = val_opt(k);
        if strcmp(OptimType,'max')
            val_opt = -val_opt;
        end
        return;
    end
end

Popt = Sselect(Popt, iv(1:Ninit));
if isfield(Sys, 'init_fun') % init_fun can modify non-uncertain parameters
    Popt = Sys.init_fun(Popt);
end
if isfield(Popt, 'init_fun')
    Popt = Popt.init_fun(Popt);
end

% max function returns the opposite of the truth value
if strcmp(OptimType,'max')
    val_opt = -val_opt;
end


end

function val = fun_max(x, Sys, phi, tspan, tau)
%% function fun_max
global Ptmp fopt traj_opt found StopWhenFound xopt

if(StopWhenFound&&~isempty(found)) %positive value found, do not need to continue
    val = -found; % optimize tries to minimize the objective function, so
    return ;          % we provide it the opposite of the truth value
end

Ptmp.pts(Ptmp.dim,1)=x;
try
    Ptmp = ComputeTraj(Sys, Ptmp, tspan);
catch  %#ok<CTCH>
    warning('SOptimProp:ComputeTraj','Error during trajectory computation when optimizing. Keep going.')
    val = inf; % do not care of the exit condition -3 for nelder-mead algo, it only happens when using global optim
    return ; % we can also set val=-inf and let the function terminates
end

val = QMITL_Eval(Sys, phi, Ptmp, Ptmp.traj(1), tau);

if(val>0)
    found = val;
end

if(val>fopt)
    fopt = val;
    traj_opt = Ptmp.traj; % we can improve that by using only Ptmp instead of traj_opt and xopt
    xopt = Ptmp.pts(Ptmp.dim,1); % as ComputeTraj launch init_fun, Ptmp.pts can be different than x
end

status = ['Robustness value: ' num2str(val) ' Current optimal: ' num2str(fopt)];
rfprintf(status);
val = -val; % optimize tries to minimize the objective function, so we
            % provide it -val instead of val
end

function val = fun_min(x, Sys, phi, tspan, tau)
%% function fun_min
global Ptmp found StopWhenFound fopt traj_opt xopt

if(StopWhenFound&&~isempty(found)) %negative value found, do not need to continue
    val = found;
    return ;
end

Ptmp.pts(Ptmp.dim,1)=x;
try
    Ptmp = ComputeTraj(Sys, Ptmp, tspan);
catch %#ok<CTCH>
    warning('SOptimProp:ComputeTraj','Error during trajectory computation when optimizing. Keep going.')
    val = inf; % do not care of the exit condition -3 for nelder-mead algo, it only happens when using global optim
    return ; % we can also let the function terminates
end

val = QMITL_Eval(Sys, phi, Ptmp, Ptmp.traj(1), tau);


if(val<0)
    found = val;
end

if(val<fopt)
    fopt = val;
    traj_opt = Ptmp.traj;
    xopt = Ptmp.pts(Ptmp.dim,1);
end

status = ['Robustness value: ' num2str(val) ' Current optimal: ' num2str(fopt)];
rfprintf(status);
end

function val = fun_zero(x, Sys, phi, tspan, tau)
%% function fun_zero
global Ptmp fopt traj_opt xopt
Ptmp.pts(Ptmp.dim,1)=x;
try
    Ptmp = ComputeTraj(Sys, Ptmp, tspan);
catch %#ok<CTCH>
    warning('SOptimProp:ComputeTraj','Error during trajectory computation when optimizing. Keep going.')
    val = inf; % do not care of the exit condition -3 for nelder-mead algo, it only happens when using global optim
    return ; % we can also let the function terminates
end

val = QMITL_Eval(Sys, phi, Ptmp, Ptmp.traj(1), tau);
status = ['Robustness value: ' num2str(val) ];
rfprintf(status);

val = abs(val);
if(val<fopt)
    fopt = val;
    traj_opt = Ptmp.traj;
    xopt = Ptmp.pts(Ptmp.dim,1);
end

end

