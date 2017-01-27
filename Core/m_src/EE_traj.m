function [X, D] = EE_traj(X0, p, k)
%
% EE_TRAJ computes a sequence of points in the parameter space for the
% Elementary Effect method
%
% Synopsis: [X, D] =  EE_traj(X0, p, k)
%
% Implements the strategy described in Global Sensitivity Analysis, A
% Primer, Saltelli et al, p113++
%
% Input:
%   -  X0 is the initial vector
%   -  p  is the number of levels (p-grid)
%   -  k  is the number of dimensions (must be equal to numel(X0) )
%
% Output:
%   - X  the dimension of X is k x (k+1)
%   - D  matrix describing wich component changes between two consecutives
%        rows in X. Dimension of D is k x k
%
% Note that X0 components should be one in
% {0, 1/(p-1), 2/(p-1), ..., 1-Delta} so that the trajectory remains in
% [0 1]^n
%
% See also SPropSensi, pRefine, EEffects

X0 = reshape(X0,1,numel(X0)); % to ensure we got a line array
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
