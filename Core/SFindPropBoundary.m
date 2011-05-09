function S = SFindPropBoundary(Sys, S, prop, tspan,tprop)
%  
%  S = SFindPropBoundary(Sys, S, prop, tspan, tprop)
%   
%  Find in S a set of points such that each points belong to a delaunay
%  simplex for which at least two vertices generate a trajectory of Sys on
%  time span tspan having a different boolean valuation of prop.
%    
%  
  
  if (~exist('tprop')) 
    tprop=0;
  end
  
  if (~isfield(Sys, 'rescale'))
    rescale = 1;
  else
    rescale = Sys.rescale;
  end
  
  S = SPurge(S);
  S = ComputeTraj(Sys,S, tspan);
  S = SEvalProp(Sys,S,prop, tprop);
  
  tri = delaunayn(rescale*S.pts(S.dim,:)');
    
  S.IsOnBoundary= zeros(1, size(S.pts,2));

  nb_prop = numel(prop);
  
  for i = 1:size(tri,1)

    % get property values for the vertices of the current simplex
    
    vertices = tri(i,:); 
    pr = S.props_values(:,vertices); 

    % Update the status of vertices with respect to the props
    traj_status = zeros(1, size(pr,2)); 
    for j=1:nb_prop 
      pri = cat(1, pr(j,:).val);
      pri = (pri(:,1)>0)';
      traj_status = traj_status+pri*2.^j;
    end

    % check if the status is uniform
    		
    status = traj_status(1);
    bool_bnd = ~isempty(find(traj_status~= status));
        
    if (bool_bnd)
      for j = 1:numel(vertices)
        S.IsOnBoundary(vertices(j))=1;         
      end
    end    
  end
  