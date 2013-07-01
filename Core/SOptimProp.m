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
%        - OptimType : 'Max' (default), 'Min' or 'Zero'.
%        - StopWhenFound : set to 1 to compute satisfaction for initial
%                          parameters in P0 then stops whenever  a positive
%                          ('Max') or negative ('Min') solution is found
%        - StopWhenFoundInit : same as above except that it does not
%                              necessarily compute all trajectories in P0
%        - Ninit   : tries the Ninit best initial pts
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
% See also SOptimPropLog Falsify
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
else
    error('SOptimProp:noMaxIter','The field opt.MaxIter is not provided.');
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

phi = QMITL_OptimizePredicates(Sys,phi); % optimization of the predicates


%% Initial values
options = optimset('MaxIter', MaxIter);

if(StopWhenFoundInit)
    nb_errors = 0; % number of ComputeTraj errors
    rfprintf_reset();
    val = zeros(1,size(P.pts,2));
    for ii = 1:size(P.pts, 2)
        Ptmp = Sselect(P,ii);
        try
            Ptmp = ComputeTraj(Sys, Ptmp, tspan);
            val(ii) = QMITL_Eval(Sys, phi, Ptmp, Ptmp.traj, tau);
        catch Me %#ok<CTCH> % in case an error occurs during computation of ComputeTraj
            warning('SOptimProp:','Error during computation of an initial trajectory, keep going.\n');
            disp(Me.getReport);
            
            if strcmp(OptimType,'max')
                val(ii) = -inf;
            else
                val(ii) = inf;
            end
            nb_errors = nb_errors + 1;
        end
        status = ['Init ' num2str(ii) '/' num2str(size(P.pts, 2)) ' Robustness value: ' num2str(val(ii))];
        rfprintf(status);
        
        switch OptimType
            case 'max'
                if(val(ii)>0)
                    val_opt = val(ii);
                    Popt = Ptmp;
                    fprintf('\n'); % to have a pretty display
                    return;
                end
            case 'min'
                if(val(ii)<0)
                    val_opt = val(ii);
                    Popt = Ptmp;
                    fprintf('\n'); % to have a pretty display
                    return;
                end
        end
        
        if(ii==1)
            Popt = Ptmp; % first time in the loop, do affectation, not concatenatation
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
    Ptmp = Sselect(Popt,1); % Initialisation of Ptmp
end


switch OptimType
    case 'max'
        [~, iv] = sort(-val);
        if(val(iv(1))>0) % if the highest value is positive
            found = val(iv(1));
        end
        fun = @(x) fun_max(x, Sys, phi, tspan, tau);
        
    case 'min'
        [~, iv] = sort(val);
        if(val(iv(1))<0) % if the lowest value is negative
            found = val(iv(1));
        end
        fun = @(x) fun_min(x, Sys, phi, tspan, tau);
        
    case 'zero'
        [~, iv] = sort(abs(val));
        fun = @(x) fun_zero(x, Sys, phi, tspan, tau);
end

if( ((StopWhenFound)&&(~isempty(found))) || (MaxIter==0) )
    Popt = Sselect(Popt,iv(1));
    val_opt = found;
    return ;
end

%% Main Loop

Ninit = min(Ninit,numel(iv));
val_opt = zeros(1,Ninit);
kk=0;
for ii = iv(1:Ninit)
    kk = kk+1;
    if isfield(opt, 'lbound')
        lbound = opt.lbound;
    else
        if ~isempty(setdiff(dim,P.dim))
            error('SOptimProp:InvalidParam',...
                'A parameter in opt.param is not in P.dim. Provide opt.lbound or modify P.');
        end
        [~,~,idim] = intersect(dim,P.dim,'stable'); % get indexes for epsi
        lbound = P.pts(dim,ii)-P.epsi(idim,ii);
    end
    
    if isfield(opt, 'ubound')
        ubound = opt.ubound;
    else
        if ~isempty(setdiff(dim,P.dim))
            error('SOptimProp:InvalidParam',...
                'A parameter in opt.param is not in P.dim. Provide opt.ubound or modify P.');
        end
        [~,~,idim] = intersect(dim,P.dim,'stable'); % get indexes for epsi
        ubound = P.pts(dim,ii)+P.epsi(idim,ii);
    end
    
    if any(lbound>ubound)
        error('SOptimProp:badIntervals','A lower bound is higher than the corresponding upper one.');
    end
    
    fprintf('\nOptimize from init point %d/%d Initial value: %g\n', kk, numel(iv), val(ii));
    rfprintf_reset();
    x0 = Popt.pts(dim,ii);
    fopt = val(ii); % we initialize with the only truth value computed for this set of values
    traj_opt = Popt.traj(Popt.traj_ref(ii));           % <--- !!! NOT SURE OF THAT (but I guess it is correct)
    xopt = Popt.pts(dim,ii);
    [~, val_opt(kk)] = optimize(fun,x0,lbound,ubound,[],[],[],[],[],[],options,'NelderMead');
    fprintf('\n');
    Popt.pts(dim,ii) = xopt;
    Popt.traj(Popt.traj_ref(ii)) = traj_opt;
    Popt.Xf(:,ii) = traj_opt.X(:,end);
    
    if((StopWhenFound)&&(~isempty(found)))
        Popt = Sselect(Popt,ii);
        if isfield(Sys, 'init_fun') % init_fun can modify non-uncertain parameters
            Popt = Sys.init_fun(Popt);
        end
        if isfield(Popt, 'init_fun')
            Popt = Popt.init_fun(Popt);
        end
        val_opt = val_opt(kk);
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

