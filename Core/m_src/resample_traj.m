function traj= resample_traj(traj, dt)
  
  if numel(dt)==1
    time= traj.time;
    new_time = time(1):dt:time(end);
  else
    time = traj.time;
      new_time = dt;
  end
  
  traj.time = new_time;
  
  newX = interp1(time,traj.X', new_time);
  
  if (size(traj.X,1)==1)
      traj.X= newX;
  else
      traj.X= newX';
  end
