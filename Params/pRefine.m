function S = pRefine(S,p,r)
% 
% PREFINE   pick initial points for the Morris global sensitivity measure
%  
% Usage:  Pr = pRefine(P, p, r)
% 
%         Refines P with r points chosen randomly in the grid with p levels
%         (See Saltelli's Books, chapter 3 or 4) Uses the suggested value
%         for Delta: p/(2(p-1)) 
%   
%         
  
% define the admissible grid and pick points

  n = numel(S.dim);           % dimension
  delta = p/(2*(p-1));        
  ngrid = floor(p*(1-delta));
  
  Stmp.dim = 1:n;
  Stmp.pts = (ngrid+1)/2*ones(n,1);
  Stmp.epsi = (ngrid+1)/2*ones(n,1);
  
  Stmp = QuasiRefine(Stmp,r);  
  Stmp.pts = floor(Stmp.pts)/(p-1);
  
  % can be improved
  S2 =Stmp; 
  S2.pts = [];
  S2.epsi = []; 
  S2.D = [];
  
% generate the elementary effects "trajectories"
  
  for k=1:r
  
    [X D] = EE_traj(Stmp.pts(:,k), p, n);
    S2.pts =  [S2.pts X];
    S2.epsi =  [S2.epsi repmat(Stmp.epsi(:,k), [1 n+1] )];
    S2.D = [S2.D D];  
    
  end
% Normalize to S ranges
    
S.pts = repmat(S.pts, [1 size(S2.pts,2)]);
S.epsi = repmat(S.epsi, [1 size(S2.pts,2)]);
S.pts(S.dim,:) = S.pts(S.dim,:) + (2*S2.pts-1).*S.epsi;
S.D = S2.D;  
 
