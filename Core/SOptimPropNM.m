function [val_opt, Popt, status, nb_call_total]  = SOptimPropNM(Sys, P, phi, opt)
%SOPTIMPROP optimizes the satisfaction of a property
% 
% Synopsis: [val_opt, Popt, status] = SOptimProp(Sys, P, phi, opt)
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
%       - MaxIter : Maximal number of optimization iterations.
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
%       - timeout : (optional), default = inf; forces convergence after 
%                   timeout s. of computation        
%       - nb_max_call : (optional), default = inf; forces convergence after 
%                   nb_max_call simulations
%
%
%    
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
%  - status  : Can be 1 (found), 2 ( stopped after max iteration), 3 (stopped after max number of simulations), -1 (stopped)
%
%  See also SOptimPropLog Falsify
%

% FIXME waitbar is cool but should be given the option to deactivate it (overhead not negligible)
   
global Ptmp; % temporary param set used to get non-variables parameter values in optim func
global found; % non empty if we found a positive or negative truth value of phi
global StopWhenFound; % See doc
global fopt; % best truth value of phi found for the current initial set of value in P
global traj_opt; % trajectory leading to fopt truth value of phi
global xopt; % parameter set leading to the trajectory traj_opt
global timeout;
global x_is_prop_param; 
global wait_bar_handle stop_requested wait_bar_step; % for the waitbar
global nb_call;
global nb_max_call; 

stop_requested = 0; 
wait_bar_handle = waitbar(0,'Optimizing formula...', 'CreateCancelBtn', 'setappdata(gcf, ''cancel'', 1)'); 

%% process options

found = [];
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
    if ~isnumeric(opt.params)
        dim = FindParam(P, opt.params);
    else
        dim = opt.params;
    end
    dim = dim(dim>0);
    dim = dim(dim<=size(P.pts,1)); % keep only existing parameters (either system or constraint parameter)
    
    P = SAddUncertainParam(P, dim);
    P = SDelUncertainParam(P, P.dim, dim);
end
dim = P.dim;

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

StopWhenFound = getfield_default(opt, 'StopWhenFound', 0);
StopWhenFoundInit = getfield_default(opt, 'StopWhenFoundInit', 0);
StopWhenFound = StopWhenFound | StopWhenFoundInit;
Ninit = getfield_default(opt, 'Ninit', size(P.pts, 2) );
timeout = getfield_default(opt, 'timeout', inf);
nb_call = 0;
nb_max_call = getfield_default(opt, 'nb_max_call', inf);

tic;
phi = STL_OptimizePredicates(Sys,phi); % optimization of the predicates

%% Initial values
if(StopWhenFoundInit)
    nb_errors = 0; % number of ComputeTraj errors
    rfprintf_reset();
    val = zeros(1,size(P.pts,2));
    for ii = 1:size(P.pts, 2)
        Ptmp = Sselect(P,ii);
        try
            Ptmp = ComputeTraj(Sys, Ptmp, tspan);
            val(ii) = STL_Eval(Sys, phi, Ptmp, Ptmp.traj, tau);
            nb_call = nb_call+1;
        catch Me % in case an error occurs during computation of ComputeTraj
            warning('SOptimProp:ComputeTraj','Error during computation of an initial trajectory, keep going.\n');
            disp(Me.getReport);
            
            if strcmp(OptimType,'max')
                val(ii) = -inf;
            else
                val(ii) = inf;
            end
            nb_errors = nb_errors + 1;
        end
        
        % Update wait bar
        wait_bar_title = sprintf('Init %d/%d Rob: %g\n',ii,size(P.pts,2),val(ii));
        set(wait_bar_handle, 'Name',wait_bar_title);
        wait_bar_status = sprintf('Current Robustness: %g',val(ii));
        wait_bar_step = ii/(size(P.pts,2)+1);
        waitbar(wait_bar_step, wait_bar_handle, wait_bar_status);
        
        drawnow; % for waitbar
        stop_requested = ~isempty(getappdata(wait_bar_handle, 'cancel'));

        if( (strcmp(OptimType,'max')&&val(ii)>0) || (strcmp(OptimType,'min')&&val(ii)<0) )
            val_opt = val(ii);
            Popt = Ptmp;
            fprintf('\n'); 
            delete(wait_bar_handle);
            status=1;
            nb_call_total= nb_call;
            return;
        end
        
        if ( stop_requested )
            val_opt = val(ii);
            Popt = Ptmp;
            fprintf('\n'); 
            delete(wait_bar_handle);
            status = -1;
            nb_call_total= nb_call;
            return;
        end

        if ( nb_call >nb_max_call) 
            val_opt = val(ii);
            Popt = Ptmp;
            fprintf('\n'); 
            delete(wait_bar_handle);
            status = 3;
            nb_call_total= nb_call;
            return;
        end
            
        if(ii==1)
            Popt = Ptmp; % first time in the loop, do affectation, not concatenatation
        else
            Popt = SConcat(Popt, Ptmp);
        end
    end
    Ninit = min(Ninit,size(P.pts,2)-nb_errors); % we dont consider parameter sets generating an error of ComputeTraj
else % FIXME timeout and nb_max_call not implemented here  
    try
        Popt = ComputeTraj(Sys, P, tspan);
    catch err
        fprintf([err.message,'\n']);
        error('SOptimProp:ComputeTraj','Error during computation of initial trajectories. Try with opt.StopWhenFoundInit=1.')
    end
    [Popt, val] = SEvalProp(Sys, Popt, phi, tau);
    Ptmp = Sselect(Popt,1); % Initialisation of Ptmp
end

% sorts initial values (most promising to least promising)
switch OptimType
    case 'max'
        [~, iv] = sort(-val);
        if(val(iv(1))>0) % if the highest value is positive
            found = val(iv(1));
        end
        fun = @(x) fun_opt(1,x, Sys, phi, tspan, tau);
        
    case 'min'
        [~, iv] = sort(val);
        if(val(iv(1))<0) % if the lowest value is negative
            found = val(iv(1));
        end
        fun = @(x) fun_opt(-1, x, Sys, phi, tspan, tau);
        
    case 'zero'
        [~, iv] = sort(abs(val));
        fun = @(x) fun_zero(x, Sys, phi, tspan, tau);
end

if ((StopWhenFound)&&(~isempty(found)))
    Popt = Sselect(Popt,iv(1));
    val_opt = val(iv(1)); 
    delete(wait_bar_handle);
    status = 1;
    nb_call_total= nb_call;
    return ;
end

if (MaxIter<=0) 
    delete(wait_bar_handle);
    status =2;
    nb_call_total= nb_call;    
    val_opt = val(iv(1)); 
    return;
end

%% Main Loop

x_is_prop_param = all(dim>Sys.DimP);

% Check param again (eventually generate idim)
if ~isfield(opt,'lbound')
    if ~isempty(setdiff(dim,P.dim))
        error('SOptimProp:InvalidParam',...
            'A parameter in opt.param is not in P.dim. Provide opt.lbound or modify P.');
    end
    [~,~,idim] = intersect(dim,P.dim,'stable'); % get indices for epsi
end
if ~isfield(opt,'ubound')
    if ~isempty(setdiff(dim,P.dim))
        error('SOptimProp:InvalidParam',...
            'A parameter in opt.param is not in P.dim. Provide opt.ubound or modify P.');
    end
    [~,~,idim] = intersect(dim,P.dim,'stable'); % get indices for epsi
end
  
options = optimset('MaxIter', MaxIter);
Ninit = min(Ninit,numel(iv));
val_opt = zeros(1,Ninit);
kk=0;
for ii = iv(1:Ninit)
    kk = kk+1;
    if isfield(opt, 'lbound')
        lbound = opt.lbound;
    else
        lbound = P.pts(dim,ii)-P.epsi(idim,ii);
    end
    
    if isfield(opt, 'ubound')
        ubound = opt.ubound;
    else
        ubound = P.pts(dim,ii)+P.epsi(idim,ii);
    end
    
    if any(lbound>ubound)
        error('SOptimProp:badIntervals','A lower bound is higher than the corresponding upper one.');
    end
    
    % update wait bar
    wait_bar_title = sprintf('Init %d/%d with rob: %g', kk, numel(iv), val(ii));
    set(wait_bar_handle, 'Name',wait_bar_title);            
    wait_bar_step = kk/(numel(iv)+1);
    
    
    x0 = Popt.pts(dim,ii);
    fopt = val(ii); % we initialize with the only truth value computed for this set of values
    traj_opt = Popt.traj{Popt.traj_ref(ii});
    xopt = Popt.pts(dim,ii);
    [~, val_opt(kk)] = optimize(fun,x0,lbound,ubound,[],[],[],[],[],[],options,'NelderMead');
    
    Popt.pts(dim,ii) = xopt;
    Popt.traj{Popt.traj_ref(ii}) = traj_opt;
    Popt.Xf(:,ii) = traj_opt.X(:,end);
    
    if((StopWhenFound)&&(~isempty(found)))
        Popt = Sselect(Popt,ii);
        if isfield(Sys, 'init_fun') % init_fun can modify non-uncertain parameters
            Popt = Sys.init_fun(Popt);
        end
        if isfield(Popt, 'init_fun')
            Popt = Popt.init_fun(Popt);
        end
        val_opt = val_opt(kk); %NM: why not found or fopt ?!
        if strcmp(OptimType,'max')
            val_opt = -val_opt;
        end
        delete(wait_bar_handle);
        status =1;
        nb_call_total= nb_call;
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

if (stop_requested)
    status = -1;
elseif (nb_call>nb_max_call)
    status=3;
else
    status=2;
end

nb_call_total= nb_call;
delete(wait_bar_handle);


end

function val = fun_opt(optim_dir, x, Sys, phi, tspan, tau)
%% function fun 
global Ptmp found StopWhenFound fopt traj_opt xopt timeout x_is_prop_param nb_call nb_max_call 
global wait_bar_handle stop_requested wait_bar_step

% updates the waitbar 
drawnow;
stop_requested = ~isempty(getappdata(wait_bar_handle, 'cancel'));

if(StopWhenFound&&~isempty(found)) %negative value found, do not need to continue
    val = found;
    return ;
end

ct = toc; 
if(ct>timeout||stop_requested)
  val = fopt;  %forces convergence of the optimizer in case timeout occured
  return;
end

nb_call = nb_call+1;
if nb_call>nb_max_call
    val = fopt;
    return;
end

Ptmp.pts(Ptmp.dim,1) = x;

if (x_is_prop_param == 0)
    try
        Ptmp = SPurge(Ptmp);
        Ptmp = ComputeTraj(Sys, Ptmp, tspan);
    catch %#ok<CTCH>
        warning('SOptimProp:ComputeTraj','Error during trajectory computation when optimizing. Keep going.')
        val = inf; % do not care of the exit condition -3 for nelder-mead algo, it only happens when using global optim
        return ; % we can also let the function terminates
    end
end

val = STL_Eval(Sys, phi, Ptmp, Ptmp.traj{1}, tau);

switch optim_dir

    case -1
        if(val<0)
            found = val;
        end

        if(val<fopt)
            fopt = val;
            traj_opt = Ptmp.traj;
            xopt = Ptmp.pts(Ptmp.dim,1);
        end
    case 1
        if(val>0)
            found = val;
        end
        
        if(val>fopt)
            fopt = val;
            traj_opt = Ptmp.traj; % we can improve that by using only Ptmp instead of traj_opt and xopt
            xopt = Ptmp.pts(Ptmp.dim,1); % as ComputeTraj launch init_fun, Ptmp.pts can be different than x
        end
        val = -val;
end

wait_bar_status = sprintf('Sim#: %g/%g  Rob: %g  Opt: %g',nb_call,nb_max_call, val,fopt);
waitbar(wait_bar_step, wait_bar_handle, wait_bar_status);

end

function val = fun_zero(x, Sys, phi, tspan, tau)
%% FIXME function fun_zero 
global Ptmp fopt traj_opt xopt timeout 

ct = toc; 
if(ct>timeout)
  val = fopt;  %forces convergence of the optimizer in case timeout occured
  return ;
end

Ptmp.pts(Ptmp.dim,1)=x;
try

    Ptmp = SPurge(Ptmp);
    Ptmp = ComputeTraj(Sys, Ptmp, tspan);
    val = STL_Eval(Sys, phi, Ptmp, Ptmp.traj{1}, tau);
catch %#ok<CTCH>
    warning('SOptimProp:ComputeTraj','Error during trajectory computation when optimizing. Keep going.')
    val = inf; % do not care of the exit condition -3 for nelder-mead algo, it only happens when using global optim
    return ; % we can also let the function terminates
end

status = sprintf('Robustness value: %g Current optimal:%g',val,fopt);
rfprintf(status);

val = abs(val);
if(val<fopt)
    fopt = val;
    traj_opt = Ptmp.traj;
    xopt = Ptmp.pts(Ptmp.dim,1);
end

end
