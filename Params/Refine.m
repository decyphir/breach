function S = Refine(S0,delta)
%      
%   REFINE Generates grid points in a N-dimensional parameter set 
%  
%   Usage: P = Refine(P0,delta)
%
%   - P0       is a parameter set
%
%   - delta    should be a scalar or a vector of dimension size(P0.dim). If it
%              is a scalar, it is interpreted as vector with all components
%              equal to its value. Its value(s) can be either integer or
%              real. If it is integer, P has delta(i) points in dimension
%              i. If it is real then delta is interpreted as a distance and
%              P is divided into points that at distance delta from one
%              another. 
%   
% 
%   This function is better understood through examples. First, create
%   some parameter set, 2 dimensions ranging from -1 to 1:
%    
%   P0 = CreateSampling(Sys, [1 2], [-1 1 -1 1]);
%
%   Then :   
%   
%   Refine(P0, 2) % creates a 2x2 grid in [-1 1 -1 1]
%    
%   Refine(P0, [2 3]  % creates a 2x3 grid 
%  
%   Refine(P0, .1) % creates a grid with resolution (.1,.1)
%  
%   Refine(P0 [.1 .2]) % creates a grid with resolution (.1,.2)
%
%   Refine(P0 [1.1 1.2]) % returns P0
%  

  n = numel(S0.dim);
  
  if (isscalar(delta))
    delta = delta*ones(numel(S0.dim),1);
    
  elseif (size(delta,1)==1)
    delta=delta';
  end
  
  if (numel(delta)>numel(S0.dim))
    delta = delta(S0.dim);
  end;
  
  S.dim = S0.dim;
  S.pts =[];
  S.epsi = [];
  
  for i = 1:size(S0.pts,2)  
    deltai=delta;
    if (sum(abs(deltai-floor(deltai))));
      deltai(deltai>2*S0.epsi(:,i))=2*S0.epsi(deltai>2*S0.epsi(:,i));
      deltai = ceil(2*S0.epsi(:,i)./deltai);
    end
  
    nb_new = prod(deltai);

    if (nb_new > 1)
    
      l = N2Nn(n,deltai);
      xlim = [ S0.pts(S0.dim,i)-S0.epsi(:,i)  S0.pts(S0.dim,i)+S0.epsi(:,i) ];
      X = repmat(S0.pts(:,i),1,nb_new);
      nepsi = repmat(S0.epsi(:,i)./(deltai),1,size(X,2));
    
      for j=1:n
        if (deltai(j)>1)
          d1 = xlim(j,1);
          d2 = xlim(j,2);
          dx(j) = (d2-d1)./(deltai(j));
          X(S0.dim(j),:) = l(j,:)*dx(j)+d1-dx(j)/2; 
        end
      end
    
      S.pts = [S.pts X];
      S.epsi = [S.epsi nepsi];
    end
  
  end

  if (isfield(S0,'traj_plot_opt'))
    S.traj_plot_opt = S0.traj_plot_opt;
  end
  
  if (isfield(S0,'X0plot_opt'))
    S.X0plot_opt = S0.X0plot_opt;
  end
  
  if (isfield(S0,'plot_proj'))
    S.plot_proj = S0.plot_proj;
  end
  
  if (isfield(S0,'ParamList'))
    S.ParamList = S0.ParamList;
  end
  
  if (isfield(S0,'selected'))
    S.selected = zeros(1, size(S.pts,2));    
  end

  S.DimX= S0.DimX;
  S.DimP= S0.DimP;


  