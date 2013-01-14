function [X, D] = EE_traj(X0, p, k)
%  
% EE_TRAJ computes a sequence of points in the parameter space for the Elementary Effect method
% 
% Synopsis: [X D] =  EE_traj(X0, p, k)
%    
%     Implements the strategy described in Global Sensitivity Analysis, A Primer, Saltelli et al, p113++
%   
%     X0 is the initial vector
%      p is the number of levels (p-grid)
%      k  is the number of dimensions
%  
%     Note that X0 components should be one in  {0, 1/(p-1), 2/(p-1), ..., 1-Delta} 
%     so that the trajectory remains in [0 1]^n

  X0 = X0';
  Delta = p/(2*(p-1)); 
  B = tril(~eye(k+1)); % triangular lower matrix
  B = B(:,1:k);
  J1 = ones(k+1,1);
  Jk = ones(k+1,k);     
  D = diag(sign(rand(1,k)-.5));
  P = diag(ones(k,1));
  P = P(:, randperm(k));
  X = (J1*X0 + (Delta/2)*((2*B-Jk)*D+Jk))*P;
  
  X = X';
  D = (D*P)';
  
end