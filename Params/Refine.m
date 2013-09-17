function P = Refine(P0, delta)
%REFINE Generates grid points in a N-dimensional parameter set
%
% Synopsis: P = Refine(P0, delta)
%
% Inputs:
%   - P0       is a parameter set
%
%   - delta    should be a scalar or a vector of dimension
%              numel(P0.dim) x 1. If it is a scalar, it is interpreted as
%              vector with all components equal to its value. If it is an
%              array, the ith line corresponds to the ith parameter in
%              P0.dim. Its value(s) can be either integer or real. If it is
%              integer, P has delta(i) points in dimension i. If it is real
%              (ie: any(floor(delta) ~= delta) is true), then delta is
%              interpreted as a distance and P is divided into points that
%              at distance delta from one another.
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
%   Refine(P0, [.1 .2]) % creates a grid with resolution (.1,.2)
%
%   Refine(P0, [2.1 2.2]) % returns P0
%
%See also CreateParamSet RandomLogRefine QuasiRefine
%


% check delta
if(delta==1)
    P = P0;
    return;
end

if isscalar(delta)
    delta = delta*ones(numel(P0.dim),1);
elseif(size(delta,1)==1)
    delta = delta'; % try to transpose if needed
end

if(numel(delta)>numel(P0.dim))
    delta = delta(P0.dim); %NM: CAUTION : may lead to error if P0.dim > numel(delta)
    %NM: better solution is
    % delta = delta(1:numel(P0.dim));
end;

n = numel(P0.dim);
P.dim = P0.dim;
P.pts = [];
P.epsi = [];

for ii = 1:size(P0.pts,2)
    deltai = delta;
    if sum(abs(deltai-floor(deltai))) % if delta contains reals
        % we compute how many points will be generated in each dimension
        deltai(deltai>2*P0.epsi(:,ii)) = 2*P0.epsi(deltai>2*P0.epsi(:,ii));
        deltai = ceil(2*P0.epsi(:,ii)./deltai);
    end
    
    nb_new = prod(deltai); % number of generated parameter sets
    
    if(nb_new > 1)
        l = N2Nn(n,deltai);
        xlim = [ P0.pts(P0.dim,ii)-P0.epsi(:,ii)  P0.pts(P0.dim,ii)+P0.epsi(:,ii) ];
        X = repmat(P0.pts(:,ii),1,nb_new);
        nepsi = repmat(P0.epsi(:,ii)./(deltai),1,size(X,2));

        for jj=1:n
            if(deltai(jj)>1)
                d1 = xlim(jj,1);
                d2 = xlim(jj,2);
                dx = (d2-d1)./(deltai(jj));
                X(P0.dim(jj),:) = l(jj,:)*dx+d1-dx/2;
            end
        end
        
        P.pts = [P.pts X]; % we cant know a priori the size of P.pts
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

P.traj_ref = IC';

if isfield(P0,'traj')
    P.traj = P0.traj;
    P.Xf = P0.Xf;
    if ~isequal(P.pts(1:P.DimP,IA), vertcat(P.traj.param))
        P = SPurge(P); % set traj_to_compute
    end
else
    P.traj_to_compute = IA';
end

end
