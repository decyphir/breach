function [val_opt Sopt val_init Sinit]  = SOptimProp(Sys, S, prop, opt)
%
% SOPTIMPROP optimizes the satisfaction of a property
%
% Synopsis: Sopt  = SOptimProp(Sys, P0, phi, tspan, lbound, ubound) 
% 
%    - S0 is a parameter set for Sys
%    - phi is a QMITL property
%    - opt is an option structure with the following fields :
%       
%        - tspan is the time domain computation of the trajectories
%        - params : variable (search) parameters    
%        - lbound : lower bounds for the search domain
%        - ubound : upper bounds for the search domain
%        - MaxIter : max number of optimization iteration 
%        - OptimType : 'Max', 'Min' or 'Zero' 
%  

%% process options

  global Stmp found
  
  found = [];
  if isfield(opt, 'tspan')
    tspan = opt.tspan;
  elseif isfield(Sys, 'tspan')
    tspan = Sys.tspan;
  elseif isfield(S, 'traj')   
    tspan = S.traj(1).time;
  else 
    tspan = 0:.2:10
  end  
  
  dim = S.dim;
        
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
  end
  
  Stmp = Sselect(S,1); 
  
  switch OptimType 
   case 'Max'
    fun = @(x) fun_max(x,Sys, Stmp, prop, tspan);       
   case 'Min'  
    fun = @(x) fun_min(x,Sys, Stmp, prop, tspan);       
   case 'Zero' 
    fun = @(x) fun_zero(x,Sys, Stmp, prop, tspan);       
  end
    
%% Initial values
  options = optimset('MaxIter', MaxIter);  
  S = ComputeTraj(Sys, S, tspan );
  [Sinit val] = SEvalProp(Sys, S, prop, 0);
  
  switch OptimType 
   case 'Max'
    [val_init iv] = sort(-val);
    Sinit = Sselect(Sinit,iv);
    fun = @(x) fun_max(x,Sys, prop, tspan);       
   
   case 'Min'  
    [val_init iv] = sort(val);
    Sinit = Sselect(Sinit,iv);
    fun = @(x) fun_min(x,Sys, prop, tspan);       
    
   case 'Zero' 
    [val_init iv] = sort(abs(val));
    Sinit = Sselect(Sinit,iv);
    fun = @(x) fun_zero(x,Sys, prop, tspan);          
  end
  
  %% Main Loop
  
  Sopt = Sinit;  
  
  for i = 1:size(S.pts, 2)
    if isfield(opt, 'lbound')
      lbound = opt.lbound;
    else
      lbound = S.pts(dim,i)-S.epsi(:,i);
    end
  
    if isfield(opt, 'ubound')
      ubound = opt.ubound;
    else
      ubound = S.pts(dim,i)+S.epsi(:,i);
    end
               
    fprintf('Optimize from init point %d/%d\n',i, size(S.pts, 2) );
    rfprintf_reset();
    x0 = S.pts(dim,i); 
    [x val_opt(i)] = optimize(fun, x0, lbound, ubound,[],[],[],[],[],[],options);
    fprintf('\n');
    Sopt.pts(dim,i) = x;    
    Sopt.traj(i) = Stmp.traj;    
  
  end
end
    
 
function val = fun_max(x, Sys, prop, tspan)
  global Stmp found
  if ~isempty(found)
    val = found;
  else
    Stmp.pts(Stmp.dim)=x;
    Stmp = ComputeTraj(Sys, Stmp, tspan);
    val = QMITL_Eval(Sys,prop, Stmp, Stmp.traj(1),0);
  end
  
  status = ['Robustness value: ' num2str(val) ];
  rfprintf(status);
  val = -val;
end
  
function val = fun_min(x, Sys, prop, tspan) 
  global Stmp found
  
  if ~isempty(found)
    val = found;
  else
    Stmp.pts(Stmp.dim)=x;
    Stmp = ComputeTraj(Sys, Stmp, tspan);
    val = QMITL_Eval(Sys,prop, Stmp, Stmp.traj(1),0);
  end
  status = ['Robustness value: ' num2str(val) ];
  rfprintf(status);
end
  
function val = fun_zero(x, Sys, prop, tspan)
  global Stmp;
  Stmp.pts(Stmp.dim)=x;
  Stmp = ComputeTraj(Sys, Stmp, tspan);
  val = QMITL_Eval(Sys,prop, Stmp, Stmp.traj(1),0);
  status = ['Robustness value: ' num2str(val) ];
  rfprintf(status);
  val = abs(val);
end
  
 
                                                      