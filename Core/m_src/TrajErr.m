function S = TrajErr(S)

%  
%  S = TrajErr(S)
%  
%  Compute the average distance between a set of trajectories (S.traj) and a set of
%  trajectories (S.etraj) estimated by sensitivity
% 
    
  nb_traj = numel(S.traj);

  N = size(S.traj{1}.X,1);
  Ns = size(S.traj{1}.XS,1)/N;
  
  for j = 1:nb_traj
    
    % Synchronize real and estimated trajectories 
    
    n = length(S.traj{j}.time);
    Xt = interp1(S.etraj(j).time', S.etraj(j).X', S.traj{j}.time')';
    XSt = interp1(S.etraj(j).time', S.etraj(j).XS', S.traj{j}.time')';    
    
    % Sensitivity matrices difference
    
    dXSt = abs(S.traj{j}.XS-XSt);

    % Construct the initial perturbation vector of the same size as the
    % trajectory  
    
    EPSI=[];
    epsi = S.epsi(:,j);
    for is = 0:Ns-1
      for i = 1:N 
        EPSI(i+N*is,1) = epsi(is+1);     
      end
    end 
    EPSI = repmat(EPSI,1,n);
    
    % Construct the expansion error vector for the trajectory
    
    EXPAErr = EPSI.*dXSt;
    ExpaErr = zeros(N,n);
    for is = 1:Ns-1
      ExpaErr = ExpaErr+EXPAErr(is*N+1:(is+1)*N,:);
    end  
    S.traj{j}.ExpaErr = ExpaErr+abs(S.traj{j}.X-Xt);
  end

  %  Err = Err/nb_traj; 
%  Err = max(Err);
  
