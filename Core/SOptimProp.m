function [val_opt Sopt]  = SOptimProp(Sys, S, prop, opt)
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

  global Stmp found StopWhenFound fopt traj_opt
  
  found = []; 
  traj_opt=[];
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
  else
     StopWhenFound = 0;
  end

  if isfield(opt,'StopWhenFoundInit')
    StopWhenFoundInit = opt.StopWhenFoundInit;
  else
    StopWhenFoundInit = opt.StopWhenFoundInit;
  end
  
  Stmp = Sselect(S,1); 
  
  switch OptimType 
   case 'Max'
    fun = @(x) fun_max(x,Sys, Stmp, prop, tspan);   
    fopt = -inf;
   case 'Min'  
    fun = @(x) fun_min(x,Sys, Stmp, prop, tspan);       
    fopt = inf;
    case 'Zero' 
    fopt  =inf;
    fun = @(x) fun_zero(x,Sys, Stmp, prop, tspan);       
  end
    
%% Initial values
  options = optimset('MaxIter', MaxIter);  
  S = ComputeTraj(Sys, S, tspan );
  [S val] = SEvalProp(Sys, S, prop, 0);
  
  switch OptimType 
   case 'Max'
    [val_init iv] = sort(-val);  
    fun = @(x) fun_max(x,Sys, prop, tspan);             
    if val(iv(1))>0
      found = val(iv(1));
    end
    
   case 'Min'  
    [val_init iv] = sort(val);
    fun = @(x) fun_min(x,Sys, prop, tspan);
    if val(iv(1))<0
      found = val(iv(1));
    end
    
   case 'Zero' 
    [val_init iv] = sort(abs(val));   
    fun = @(x) fun_zero(x,Sys, prop, tspan);          
  end
  
  %% Main Loop
  
  Sopt = S;  
  k=0;
  for i = iv
    k = k+1;
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
               
    fprintf('Optimize from init point %d/%d Initial value: %g\n',k, numel(iv), val(i) );
    rfprintf_reset();
    x0 = S.pts(dim,i); 
    [x val_opt(k)] = optimize(fun, x0, lbound, ubound,[],[],[],[],[],[],options);
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
  global fopt traj_opt found StopWhenFound
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
  
 
                                                      