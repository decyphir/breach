function [val_opt, Sopt]  = SOptimProp(Sys, P, prop, opt)
%
% SOPTIMPROP optimizes the satisfaction of a property
%
% Synopsis: [val_opt, Sopt]  = SOptimProp(Sys, P0, phi, opt)
%
%    - P0 is a parameter set for Sys
%    - phi is a QMITL property
%    - opt is an option structure with the following fields :
%
%        - tspan is the time domain computation of the trajectories
%        - params : (mandatory) variable (search) parameters
%        - lbound : (mandatory) lower bounds for the search domain
%        - ubound : (mandatory) upper bounds for the search domain
%        - MaxIter : (mandatory) max number of optimization iteration
%        - OptimType : 'Max' (default), 'Min' or 'Zero'
%        - StopWhenFound : set to 1 to compute satisfaction for initial parameters
%                          in P0 then stops whenever  a positive ('Max') or
%                          negative ('Min') solution is found
%        - StopWhenFoundInit : same as above except that it does not necessarily
%                              compute all trajectories in P0
%

%% process options

global Stmp found StopWhenFound fopt traj_opt

found = [];
traj_opt=[];
if isfield(opt, 'tspan')
    tspan = opt.tspan;
elseif isfield(Sys, 'tspan')
    tspan = Sys.tspan;
elseif isfield(P, 'traj')
    tspan = P.traj(1).time;
else
    tspan = 0:.2:10;
end

dim = P.dim;

if isfield(opt,'OptimType')
    OptimType = opt.OptimType;
else
    OptimType = 'Max';
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
else
    StopWhenFoundInit = 0;
end

Stmp = Sselect(P,1);

switch OptimType
    case 'Max'
        fun = @(x) fun_max(x, Sys, Stmp, prop, tspan);
        fopt = -inf;
    case 'Min'
        fun = @(x) fun_min(x, Sys, Stmp, prop, tspan);
        fopt = inf;
    case 'Zero'
        fopt = inf;
        fun = @(x) fun_zero(x, Sys, Stmp, prop, tspan);
end

%% Initial values
options = optimset('MaxIter', MaxIter);

if (StopWhenFoundInit)
    rfprintf_reset();
    val = zeros(1,size(P.pts,2));
    for i = 1:size(P.pts, 2)
        Stmp = Sselect(P,i);
        Stmp = ComputeTraj(Sys, Stmp, tspan );
        val(i) = QMITL_Eval(Sys, prop, Stmp, Stmp.traj, 0);
        status = ['Init ' num2str(i) '/' num2str(size(P.pts, 2)) ' Robustness value: ' num2str(val(i)) ];
        rfprintf(status);
        
        switch OptimType
            case 'Max'
                if val(i)>0
                    val_opt = val(i);
                    Sopt = Stmp;
                    return;
                end
                
            case 'Min'
                if val(i)<0
                    val_opt = val(i);
                    Sopt = Stmp;
                    return;
                end
        end
        
        if i==1
            Sopt = Stmp;
        else
            Sopt = SConcat(Sopt, Stmp);
        end
    end
else
    Sopt = ComputeTraj(Sys, P, tspan);
    [Sopt, val] = SEvalProp(Sys, Sopt, prop, 0);
end


switch OptimType
    case 'Max'
        [val_init, iv] = sort(-val);
        fun = @(x) fun_max(x, Sys, prop, tspan);
        if val(iv(1))>0
            found = val(iv(1));
        end
        
    case 'Min'
        [val_init, iv] = sort(val);
        fun = @(x) fun_min(x, Sys, prop, tspan);
        if val(iv(1))<0
            found = val(iv(1));
        end
        
    case 'Zero'
        [val_init, iv] = sort(abs(val));
        fun = @(x) fun_zero(x, Sys, prop, tspan);
end

%% Main Loop
val_opt = val(iv(1));

if (MaxIter==0)
    return;
end

k=0;
for i = iv
    k = k+1;
    if isfield(opt, 'lbound')
        lbound = opt.lbound;
    else
        lbound = P.pts(dim,i)-P.epsi(:,i);
    end
    
    if isfield(opt, 'ubound')
        ubound = opt.ubound;
    else
        ubound = P.pts(dim,i)+P.epsi(:,i);
    end
    
    fprintf('\nOptimize from init point %d/%d Initial value: %g\n',k, numel(iv), val(i) );
    rfprintf_reset();
    x0 = P.pts(dim,i);
    [x, val_opt(k)] = optimize(fun, x0, lbound, ubound,[],[],[],[],[],[],options);
    fprintf('\n');
    Sopt.pts(dim,i) = x;
    Sopt.traj(Sopt.traj_ref(i)) = traj_opt;
    Sopt.Xf(:,i) = traj_opt.X(:,end);
    
    if (StopWhenFoundInit)&&(~isempty(found))
        Sopt = Sselect(Sopt,i);
        val_opt = val_opt(k);
        break
    end
end




end

function val = fun_max(x, Sys, prop, tspan)
global Stmp fopt traj_opt found StopWhenFound
if (~isempty(found)&&StopWhenFound)
    val = found;
else
    Stmp.pts(Stmp.dim)=x;
    Stmp = ComputeTraj(Sys, Stmp, tspan);
    val = QMITL_Eval(Sys,prop, Stmp, Stmp.traj(1),0);
end

if (val>0)
    found = val;
end

if (val>fopt)
    fopt = val;
    traj_opt = Stmp.traj;
end

status = ['Robustness value: ' num2str(val) 'Current optimal: ' num2str(fopt)];
rfprintf(status);
val = -val;
end

function val = fun_min(x, Sys, prop, tspan)
global Stmp found StopWhenFound fopt traj_opt

if (~isempty(found)&&StopWhenFound)
    val = found;
else
    Stmp.pts(Stmp.dim)=x;
    Stmp = ComputeTraj(Sys, Stmp, tspan);
    val = QMITL_Eval(Sys,prop, Stmp, Stmp.traj(1),0);
end

if (val<0)
    found = val;
end

if (val<fopt)
    fopt = val;
    traj_opt = Stmp.traj;
end

status = ['Robustness value: ' num2str(val) '  Current optimal: ' num2str(fopt)];
rfprintf(status);
end

function val = fun_zero(x, Sys, prop, tspan)
global Stmp fopt traj_opt
Stmp.pts(Stmp.dim)=x;
Stmp = ComputeTraj(Sys, Stmp, tspan);
val = QMITL_Eval(Sys,prop, Stmp, Stmp.traj(1),0);
status = ['Robustness value: ' num2str(val) ];
rfprintf(status);
if (abs(val)<fopt)
    fopt = abs(val);
    traj_opt = Stmp.traj;
end

val = abs(val);


end


