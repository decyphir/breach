function VoxelTraj(traj,c,alpha,dt)
  
  t = 0;
  k = 1;

  hold on;

  while(1)
      
    DX = traj.Expa(1:3,k)';
    X = traj.X(1:3,k)'-DX;
    voxel(X,2*DX,c,alpha);

    if t>traj.time(end)
      break
    end

    t = t+dt;
    k=find(traj.time>=t,1);

    if isempty(k)
      k = numel(traj.time);
    end
    
  end
    
    
