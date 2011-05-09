function Sopt = SOptimProp(Sys, S, prop, tspan, lbound, ubound)
%
% SOPTIMPROP optimizes the satisfaction of a property
%
% Usage: Sopt  = SOptimProp(Sys, S0, phi, tspan) 
% 
%    - S0 is a parameter set for Sys
%    - phi is obviously a QMITL property
%    - tspan is the time domain computation of the trajectories
%

  
  dim = S.dim;
  Stmp = select(S,1); 
  fun = @(x) fun0(x,Sys, Stmp, prop, tspan);
   
  Sopt = S;
  
  for i = 1:size(S.pts, 2)
    %    ubound = S.pts(dim,i)+S.epsi(:,i);
    %    lbound = S.pts(dim,i)-S.epsi(:,i);

    
    disp('\n\n ---------------------- New point ------------------- \n\n')
    x0 = S.pts(dim,i); 
    x = optimize(fun, x0, lbound, ubound);
    Sopt.pts(dim,i) = x;
  
    
  end
    
    
function val = fun0(x, Sys,Stmp, prop, tspan)
  
  Stmp.pts(Stmp.dim)=x;
  Sf = ComputeTraj(Sys, Stmp, tspan);
  val = QMITL_Eval(Sys,prop, Sf.traj(1),0);
%  val = abs(val);
val = -val;

%  fprintf('x(1): %g x(2): %g val: %g   \n', x(1), x(2), -val);
  fprintf('x(1): %g x(2): %g x(3): %g val: %g   \n', x(1), x(2), x(3), -val);
  
                                                      