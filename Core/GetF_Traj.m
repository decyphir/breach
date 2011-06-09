function f = GetF_Traj(Pf,ind_traj, t)

  global Sys;
  InitSystem(Sys);
  
  traj = Pf.traj(ind_traj);  

  if (exist('t'))
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