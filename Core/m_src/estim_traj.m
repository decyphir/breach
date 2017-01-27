function  etraj = estim_traj(traj, epsi)
%  
% etraj = estim_traj(traj, epsi)
%
% Computes an estimated trajectory using sensitivity matrice in traj
%  
        
  etraj.time = traj.time;
  l = numel(traj.time);
  
  N = size(traj.X,1);
  Ns = size(traj.XS,1)/N;
  
  for is = 0:Ns-1
    for i = 1:N 
      EPSI(i+N*is,1) = epsi(is+1);     
    end
  end 

  EPSI = repmat(EPSI,1,l);
  DX = EPSI.*traj.XS;
  dX = DX(1:N, :);
  
  for is = 1:Ns-1
    dX = dX+DX(is*N+1:(is+1)*N,:);
  end  
 
  etraj.X = traj.X+dX;
  etraj.XS = traj.XS;
