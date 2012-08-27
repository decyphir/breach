function Sopt = SOptimProp(Sys, S, prop, tspan, lbound, ubound)
%
% SOPTIMPROP optimizes the satisfaction of a property
%
% Synopsis: Sopt  = SOptimProp(Sys, P0, phi, tspan, lbound, ubound) 
% 
%    - S0 is a parameter set for Sys
%    - phi is a QMITL property
%    - tspan is the time domain computation of the trajectories
%    - lbound : lower bounds for the search domain
%    - ubound : upper bounds for the search domain
%
  
  dim = S.dim;
  Stmp = Sselect(S,1); 
  fun = @(x) fun0(x,Sys, Stmp, prop, tspan);
   
  Sopt = S;
  
  options = optimset('MaxIter', 100);
  for i = 1:size(S.pts, 2)
    % ubound = S.pts(dim,i)+S.epsi(:,i);
    % lbound = S.pts(dim,i)-S.epsi(:,i);
    
    disp('\n\n ---------------------- New point ------------------- \n\n')
    x0 = S.pts(dim,i); 
    x = optimize(fun, x0, lbound, ubound,[],[],[],[],[],[],options);
    Sopt.pts(dim,i) = x;
  
  end
    
    
function val = fun0(x, Sys,Stmp, prop, tspan)
  
  Stmp.pts(Stmp.dim)=x;
  Sf = ComputeTraj(Sys, Stmp, tspan);
  val = QMITL_Eval(Sys,prop, Sf.traj(1),0);
% val = abs(val);
  val = -val;

% fprintf('x(1): %g x(2): %g val: %g   \n', x(1), x(2), -val);
  fprintf('x: ');
  for i= 1:numel(x)
    fprintf('%g  ', x(i));
  end
  fprintf(' val: %g \n', -val);
  
                                                      