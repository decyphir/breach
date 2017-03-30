function P = RefineEstim(P0, delta)
%REFINEESTIM return a sampling P refining P0 depending on the parameter
% delta and compute estimated trajectories for each point in the
% refinement.
%
% Synopsis: P = RefineEstim(P0, delta)
% 
% Inputs:
%   delta should be a scalar or a vector of dimension size(P0.dim). It
%   can be either integer or real. If it is integer, P  has delta(i)
%   points in dimension i. If it is real, then P has dispersion
%   delta. If P0 had already a dispersion less than delta, it is
%   unchanged.
%


if ~isfield(P0,'traj')
    error('RefineEstim:noTrajField','Compute trajectories first.')
elseif ~isfield(P0.traj{1},'XS')
    error('RefineEstim:noXSField','Compute sensitivities first.')
end

n = numel(P0.dim);
if isscalar(delta)
    delta = delta*ones(numel(P0.dim),1);
elseif(size(delta,1)==1)
    delta = delta';
end

if(numel(delta)>numel(P0.dim))
    delta = delta(P0.dim);
end;

P.dim = P0.dim;
P.pts = [];
P.epsi = [];
P.etraj = [];

for ii = 1:size(P0.pts,2)
    deltai = delta;
    if(sum(abs(deltai-floor(deltai))));
        deltai(deltai>2*P0.epsi(:,ii)) = 2*P0.epsi(deltai>2*P0.epsi(:,ii));
        deltai = ceil(2*P0.epsi(:,ii)./deltai);
    end
    
    nb_new = prod(deltai);
    
    if(nb_new > 1)
        
        l = N2Nn(n,deltai);
        xlim = [ P0.pts(P0.dim,ii)-P0.epsi(:,ii)  P0.pts(P0.dim,ii)+P0.epsi(:,ii) ];
        X = repmat(P0.pts(:,ii),1,nb_new);
        nepsi = repmat(P0.epsi(:,ii)./(deltai),1,size(X,2));
        
        for jj=1:n
            if(deltai(jj)>1)
                d1 = xlim(jj,1);
                d2 = xlim(jj,2);
                dx(jj) = (d2-d1)./(deltai(jj));
                X(P0.dim(jj),:) = l(jj,:)*dx(jj)+d1-dx(jj)/2;
            end
        end
        P.pts = [P.pts X];
        P.epsi = [P.epsi nepsi];
        
%        etraj = zeros(1,nb_new); % initialise etraj
        for jj = 1:nb_new
            dx = X(:,jj)-P0.pts(:,ii);
            etraj(jj) = estim_traj(P0.traj{ii}, dx(P0.dim));
        end
        
        P.etraj = [P.etraj etraj];
    end
end


%   try
%       rmfield(P,'XS0');
%   end
%
if isfield(P0,'traj_plot_opt')
    P.traj_plot_opt = P0.traj_plot_opt;
end

if isfield(P0,'X0plot_opt')
    P.X0plot_opt = P0.X0plot_opt;
end

if isfield(P0,'plot_proj')
    P.plot_proj = P0.plot_proj;
end

P.DimX = P0.DimX;
P.DimP = P0.DimP;
P.ParamList = P0.ParamList;

end
