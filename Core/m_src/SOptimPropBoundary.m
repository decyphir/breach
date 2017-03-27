function S = SOptimPropBoundary(Sys, S, prop, tspan,optim_opt)
  
  
  dim = S.dim;
  Stmp = select(S,1); 
  fun = @(x) fun0(x,Sys, Stmp, prop, tspan);
  
  for i = 1:size(S.pts, 2)
    ubounds = S.pts(dim,i)+S.epsi(:,i);
    lbound = S.pts(dim,i)-S.epsi(:,i);
    x0 = S.pts(dim,i);
%    x = optimize(fun, x0, lbound, ubounds,[],[],[],[],[],...
%         [], optim_opt);
    x = optimize(fun, x0,[] ,[],[],[],[],[],[],...
         [], optim_opt);
    S.pts(dim,i) = x;
  end
    
    
function val = fun0(x, Sys,Stmp, prop, tspan)
  
  Stmp.pts(Stmp.dim)=x;
  Sf = ComputeTraj(Sys, Stmp, tspan);
  val = STL_Eval(Sys, prop,Sf.traj{1},0);
  val = -val;
  
  
