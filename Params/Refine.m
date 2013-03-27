function P = Refine(P0, delta)
%REFINE Generates grid points in a N-dimensional parameter set
%
% Synopsis: P = Refine(P0, delta)
%
% Inputs:
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
%See also CreateParamSet RandomLogRefine
%

if(delta==1)
    P = P0;
    return;
end

n = numel(P0.dim);

if isscalar(delta)
    delta = delta*ones(numel(P0.dim),1);
    
elseif(size(delta,1)==1)
    delta=delta';
end

if(numel(delta)>numel(P0.dim))
    delta = delta(P0.dim);
end;

P.dim = P0.dim;
P.pts =[];
P.epsi = [];

for i = 1:size(P0.pts,2)
    deltai=delta;
    if (sum(abs(deltai-floor(deltai))));
        deltai(deltai>2*P0.epsi(:,i)) = 2*P0.epsi(deltai>2*P0.epsi(:,i));
        deltai = ceil(2*P0.epsi(:,i)./deltai);
    end
    
    nb_new = prod(deltai);
    
    if(nb_new > 1)
        
        l = N2Nn(n,deltai);
        xlim = [ P0.pts(P0.dim,i)-P0.epsi(:,i)  P0.pts(P0.dim,i)+P0.epsi(:,i) ];
        X = repmat(P0.pts(:,i),1,nb_new);
        nepsi = repmat(P0.epsi(:,i)./(deltai),1,size(X,2));
        
        for j=1:n
            if(deltai(j)>1)
                d1 = xlim(j,1);
                d2 = xlim(j,2);
                dx(j) = (d2-d1)./(deltai(j));
                X(P0.dim(j),:) = l(j,:)*dx(j)+d1-dx(j)/2;
            end
        end
        
        P.pts = [P.pts X];
        P.epsi = [P.epsi nepsi];
    end
    
end

if isfield(P0,'traj_plot_opt')
    P.traj_plot_opt = P0.traj_plot_opt;
end

if isfield(P0,'X0plot_opt')
    P.X0plot_opt = P0.X0plot_opt;
end

if isfield(P0,'plot_proj')
    P.plot_proj = P0.plot_proj;
end

if isfield(P0,'ParamList')
    P.ParamList = P0.ParamList;
end

if isfield(P0,'selected')
    P.selected = zeros(1, size(P.pts,2));
end

P.DimX = P0.DimX;
P.DimP = P0.DimP;

%  Checks for pts sharing the sames systems parameters

X = P.pts(1:P.DimP,:)';
[~,IA,IC] = unique(X,'rows');

P.traj_ref= IC';

if (isfield(P0,'traj'))
    
    P.traj = P0.traj;
    P.Xf = P0.Xf;
    if ~isequal(P.pts(1:P.DimP,IA), vertcat(P.traj.param))
        P.traj_to_compute = IA';
    end
else
    P.traj_to_compute = IA';
end

end
