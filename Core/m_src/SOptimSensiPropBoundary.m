function S = SOptimSensiPropBoundary(Sys, S, prop, tspan,optim_opt)
  
  dim = S.dim;
  Stmp = select(S,1); 
  fun = @(x) fun0(x,Sys, Stmp, prop, tspan);
  
  for i = 1:size(S.pts, 2)

    x0 = S.pts(dim,i);
    x = minimize(x0,fun);
    S.pts(dim,i) = x;
  end
    
    
function [val vald] = fun0(x, Sys,Stmp, prop, tspan)
  
  Stmp.pts(Stmp.dim)=x;
  Sf = ComputeTrajSensi(Sys, Stmp, tspan);
  [val vald] = STL_EvalSensi(prop,Sf.traj{1},0);
  val = abs(val);
  vald= sign(val)*vald;
  
