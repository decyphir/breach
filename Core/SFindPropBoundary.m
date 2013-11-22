function Pf = SFindPropBoundary(Sys, P, phis, tspan, tprop, iter)
%SFINDPROPBOUNDARY finds parameters in a set which are close to a boundary
% between different satisfaction regions
%
%  Synopsis: Pf = SFindPropBoundary(Sys, P, prop, tspan, tprop, iter)
%
%  Find in P a set of points such that each points belong to a delaunay
%  simplex for which at least two vertices generate a trajectory of Sys on
%  time span tspan having a different boolean valuation of prop.
%
%

needs_sensi = check_sensi(phis);

if ~exist('tprop','var')
    tprop = 0;
end

if ~exist('iter','var')
    iter = 1; % will only estimate (select) the boundary
end

if(iter==0)
    return;
end

%P = SPurge(P);
%P = SPurge_props(P);

Pf = SCopyEmpty(P);
Pf.selected = [];

Pf.traj = [];
Pf.props_values = [];

nb_phis = numel(phis);

for kk=1:iter
    if needs_sensi
        P = ComputeTrajSensi(Sys, P, tspan);
    else
        P = ComputeTraj(Sys, P, tspan);
    end
    
    P = SEvalProp(Sys, P, phis, tprop);
    
    tri = delaunayn(P.pts(P.dim,:)', {'Qt', 'Qbb','Qc', 'Qz'});
    P.IsOnBoundary = zeros(1, size(P.pts,2));
    
    for ii = 1:size(tri,1)
        
        % get property values for the vertices of the current simplex
        
        vertices = tri(ii,:);
        pr = P.props_values(:,vertices);
        
        % Update the status of vertices with respect to the props
        
        traj_status = zeros(1, size(pr,2));
        for jj=1:nb_phis
            pri = cat(1, pr(jj,:).val);
            pri = (pri(:,1)>0)';
            traj_status = traj_status+pri*2.^jj;
        end
        
        % check if the status is uniform
        
        status = traj_status(1);
        bool_bnd = ~isempty(find(traj_status~=status,1));
        
        if(bool_bnd)
            for jj = 1:numel(vertices)
                P.IsOnBoundary(vertices(jj)) = 1;
            end
        end
    end
    
    Pn = Sselect(P, find(P.IsOnBoundary>0));
    idx_nn = find(P.IsOnBoundary==0);
    Pnn = Sselect(P, idx_nn);
    
    nbnn = numel(idx_nn);
    if(nbnn)
        Pf = SConcat(Pf, Pnn);
        Pf.selected(end+1:end+nbnn) = 0 ;
    end
    
    if(kk < iter)
        nbn = size(Pn.pts,2);
        if(nbn>0)
            %P = VoronoiRefine(Pn);
            P = Refine(Pn,2);
        else
            break
        end
    end
end

nbn = size(Pn.pts,2);

if(nbn > 0)
    Pf = SConcat(Pf, Pn);
    Pf.selected(end+1:end+nbn) = 1;
else
    return;
end

end
