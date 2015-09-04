function traj = s2traj(s1,s2)
  
  if (s1.t(1) < s2.t(1))
    s2.t = [s1.t(1) s2.t];
    s2.f = [s2.f(1) s2.f];
    s2.a = [0 s2.a];     
  elseif (s1.t(1) > s2.t(1))
    s1.t = [s2.t(1) s1.t];
    s1.f = [s1.f(1) s1.f];
    s1.a = [0 s1.a];
  end

  if (s1.t(end) > s2.t(end))
    s2.t = [s2.t s1.t(end) ];
    s2.f = [s2.f s2.f(end)];
    s2.a = [s2.a 0];     
  elseif (s1.t(end) < s2.t(end))
    s1.t = [s1.t s2.t(end) ];
    s1.f = [s1.f s1.f(end) ];
    s1.a = [s1.a 0];
  end
  
  dt=0.001;  
  traj.time= s1.t(1):dt:s1.t(end);
  traj.X(1,:) = interp1(s1.t, s1.f, traj.time);
  traj.X(2,:) = interp1(s2.t, s2.f, traj.time);
  
