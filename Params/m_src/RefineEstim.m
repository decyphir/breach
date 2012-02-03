function S = RefineEstim(S0,delta)

%
%   S = RefineEstim(S0,delta)
% 
%   Return a sampling S refining S0 depending on the parameter delta
%   and compute estimated trajectories for each point in the refinement.
% 
%   delta should be a scalar or a vector of dimension size(S0.dim). It
%   can be either integer or real. If it is integer, S  has delta(i)
%   points in dimension i. If it is real, then S has dispersion
%   delta. If S0 had already a dispersion less than delta, it is
%   unchanged.  
%  
  
    
  if (~isfield(S0,'traj'))
    error('Compute trajectories first')      
  elseif (~isfield(S0.traj(1),'XS'))
    error('Compute Sensitivities first')
  end      
       
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
  S.etraj = [];
    
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
      
      for j = 1:nb_new
        dx = X(:,j)-S0.pts(:,i);
        etraj(j) = estim_traj(S0.traj(i), dx(S0.dim));
      end
      
      S.etraj = [S.etraj etraj];
    end      
  end
 
    
%   try
%       rmfield(S,'XS0');    
%   end
%   
  if (isfield(S0,'traj_plot_opt'))
    S.traj_plot_opt = S0.traj_plot_opt;
  end
  
  if (isfield(S0,'X0plot_opt'))
    S.X0plot_opt = S0.X0plot_opt;
  end
  
  if (isfield(S0,'plot_proj'))
    S.plot_proj = S0.plot_proj;
  end
  
  S.DimX = S0.DimX;
  S.DimP = S0.DimP;
  S.ParamList = S0.ParamList;