function J = GetJf_Traj(Sys,Pf,ind_traj, t)

  InitSystem(Sys);
  
  traj = Pf.traj{ind_traj};  
  ind_t= find(traj.time>= t,1);
  
  t0 = traj.time(ind_t-1);
  t1 = traj.time(ind_t);
  X = interp1([t0 t1], [traj.X(:,ind_t) traj.X(:, ind_t+1)]', t);
  
  pts = Pf.pts(:, ind_traj);
  pts(1:numel(X)) = X';
  
  J = cvm(58,t,pts);
