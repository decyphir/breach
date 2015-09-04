function f = GetF_Traj(Sys, traj , t)
% GETF_TRAJ returns the values of the rhs ODE function along a computed trajectory 
%  
% Synopsys:    f = GetF_Traj(Sys, traj [, t])
%
  
  InitSystem(Sys);  
 
  if (exist('t','var'))
    ind_t= find(traj.time>= t,1);
  
    t0 = traj.time(ind_t-1);
    t1 = traj.time(ind_t);
    X = interp1([t0 t1], [traj.X(:,ind_t-1) traj.X(:, ind_t)]', t);
    
    pts = Pf.pts(:, ind_traj);
    pts(1:numel(X)) = X';
  
  else
    t = 0;
    nb_stp= numel(traj.time);
    pts = repmat(Pf.pts(:, ind_traj),1,nb_stp) ;
    pts(1:Sys.DimX, :) = traj.X;
  end
  
  f = cvm(59,t,pts);
