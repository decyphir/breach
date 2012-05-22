function S = SFindPropBoundaReach(Sys, S, prop, tspan)
%  
%  S = SFindPropBoundaReach(Sys, S, prop, tspan, tprop)
%   
%  Find in S a set of points such that the local reachable set from these
%  points intersects with the satisfaction boundary of the property prop
%  on time tspan
%    
%  
  
  if (~exist('tprop')) 
    tprop=0;
  end
        
  S = SPurge(S);
  S = ComputeTrajSensi(Sys,S, tspan);
  
  S = SEvalProp(Sys,S,prop, tprop);
  
  St = RefineEstim(S,2);
  St.traj = St.etraj; 
  
  Sp = SEvalProp(Sys,St,prop, tprop);
  mult = 2^(numel(S.dim));
  
  S.IsOnBoundary= zeros(1, size(S.pts,2));
  
  for i = 0:numel(S.traj)-1
    vertices = mult*i+1:mult*(i+1);
    bool_vec = cat(1,Sp.props_values(1,vertices).val)';
    bool_bnd = (~all(bool_vec>0))&&(any(bool_vec>0));
    if (bool_bnd)
      for j = 1:numel(vertices)
        S.IsOnBoundary(vertices(j))=1;         
      end
    end    
  end
  
  
  
  
  