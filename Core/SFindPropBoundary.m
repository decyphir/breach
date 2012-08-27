function Sf = SFindPropBoundary(Sys, S, prop, tspan,tprop, iter)
% SFINDPROPBOUNDARY  find parameters in a set which are close to a boundary
% between different satisfaction regions 
%
%  Synopsis: 
%              S = SFindPropBoundary(Sys, S, prop, tspan, tprop, iter)
%   
%  Find in S a set of points such that each points belong to a delaunay
%  simplex for which at least two vertices generate a trajectory of Sys on
%  time span tspan having a different boolean valuation of prop.
%    
%  
  
  needs_sensi = check_sensi(prop);
  
  if (~exist('tprop')) 
    tprop=0;
  end
  
  if (~exist('iter'))
    iter = 1; % will only estimate (select) the boundary
  end

  if iter ==0
    return;
  end
  
  %S = SPurge(S);
  %S = SPurge_props(S);
  
  Sf = SCopyEmpty(S);
  Sf.selected = [];

  Sf.traj = [];
  Sf.props_values = [];
  
  nb_prop = numel(prop);
  
  for k=1:iter
    if needs_sensi 
      S = ComputeTrajSensi(Sys,S, tspan);
    else
      S = ComputeTraj(Sys,S, tspan);      
    end
      
    S = SEvalProp(Sys,S,prop, tprop);
  
    tri = delaunayn(S.pts(S.dim,:)', {'Qt', 'Qbb','Qc', 'Qz'});
    S.IsOnBoundary= zeros(1, size(S.pts,2));
    
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
  
    Sn =  Sselect(S, find(S.IsOnBoundary>0));
    Snn = Sselect(S, find(S.IsOnBoundary==0));
  
    nbnn = size(Snn.pts,2);   
    if (nbnn)
      Sf = SConcat(Sf, Snn);
      Sf.selected(end+1:end+nbnn) = 0 ;   
    end

    if k < iter
      nbn = size(Sn.pts,2);
      if (nbn>0)
        %S = VoronoiRefine(Sn);
        S = Refine(Sn,2);
      else 
        break
      end
    end
  end
  
  nbn = size(Sn.pts,2);   

  if (nbn > 0)
    Sf = SConcat(Sf, Sn);
    Sf.selected(end+1:end+nbn) = 1;
  else
    return;
  end