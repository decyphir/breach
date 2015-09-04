function P = SFindPropBoundaReach(Sys, P, prop, tspan)
%  
%  P = SFindPropBoundaReach(Sys, P, prop, tspan, tprop)
%   
%  Find in P a set of points such that the local reachable set from these
%  points intersects with the satisfaction boundary of the property prop
%  on time tspan
%    
%  

if ~exist('tprop','var')
    tprop = 0;
end

P = SPurge(P);
P = ComputeTrajSensi(Sys,P, tspan);

P = SEvalProp(Sys, P, prop, tprop);

Pt = RefineEstim(P,2);
Pt.traj = Pt.etraj;

Sp = SEvalProp(Sys,Pt,prop, tprop);
mult = 2^(numel(P.dim));

P.IsOnBoundary = zeros(1, size(P.pts,2));

for ii = 0:numel(P.traj)-1
    vertices = mult*ii+1:mult*(ii+1);
    bool_vec = cat(1,Sp.props_values(1,vertices).val)';
    bool_bnd = (~all(bool_vec>0))&&(any(bool_vec>0));
    if(bool_bnd)
        for jj = 1:numel(vertices)
            P.IsOnBoundary(vertices(jj))=1;
        end
    end
end

end  
  
