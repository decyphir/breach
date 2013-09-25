function P = Refine(P0, delta)
%REFINE Generates grid points in a N-dimensional parameter set
%
% Synopsis: P = Refine(P0, delta)
%
% Inputs:
%   - P0       is a parameter set which may contain many parameter vectors
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
% Output:
%  - P  is the generated parameter set
%
% This function is better understood through examples. First, create
% some parameter set, 2 dimensions ranging from -1 to 1:
%
% Examples (lorentz84):
%   CreateSystem;
%   P0 = CreateParamSet(Sys, [1, 2], [0, 2; 0, 2]);
%
%   %Then :
%
%   P = Refine(P0, 2); % creates a 2x2 grid in [0 2 0 2]
%   P.pts   % show the generated parameter vectors
%
%   P = ComputeTraj(Sys,P0,1:0.1:10);
%   P = Refine(P, [3 3]);  % creates a 3x3 grid
%   P.traj_ref     % should be [0 0 0 0 1 0 0 0 0]
%   P.traj_to_compute  % should be [1 2 3 4 6 7 8 9]
%   P.traj
%
%   P = Refine(P0, .1); % creates a grid with resolution (.1,.1)
%   P.pts(:,[1,2,3,19,20,21,22,400])
%
%   Refine(P0, [.1 .2]) % creates a grid with resolution (.1,.2) (works
%                       % well as .1 and .2 are divisor of 2*P.epsi)
%
%   P = Refine(P0, [1.3,1.3]); % creates a grid with "resolution (1.3,1.3)"
%   P.pts    % pay attention here to the values !
%
%   P = Refine(P0, [2.1 2.2]);   % returns P0:
%   all(size(P.pts)==size(P0.pts)) & all((P.pts==P0.pts)  % true
%
%See also CreateParamSet RandomLogRefine QuasiRefine SAddUncertainParam
%SetEpsi
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
    delta = delta(P0.dim);
end;

n = numel(P0.dim);
P.dim = P0.dim;
P.pts = [];
P.epsi = [];

for ii = 1:size(P0.pts,2)
    deltai = delta;
    if sum(abs(deltai-floor(deltai))) % if delta contains reals
        % we compute how many points will be generated in each dimension
        deltai(deltai>2*P0.epsi(:,ii)) = 2*P0.epsi(deltai>2*P0.epsi(:,ii),ii);
        deltai = ceil(2*P0.epsi(:,ii)./deltai);
    end
    
    nb_new = prod(deltai); % number of generated parameter vector
    
    if(nb_new > 1)
        l = N2Nn(n, deltai);
        xlim = [ P0.pts(P0.dim,ii)-P0.epsi(:,ii) , P0.pts(P0.dim,ii)+P0.epsi(:,ii) ];
        X = repmat(P0.pts(:,ii), 1, nb_new);
        nepsi = repmat(P0.epsi(:,ii)./(deltai), 1, nb_new);

        for jj=1:n
            if(deltai(jj)>1) %NM: there maybe a bug here (delta(P0.dim(jj)) ?
                d1 = xlim(jj,1);
                d2 = xlim(jj,2);
                dx = (d2-d1)./(deltai(jj));
                X(P0.dim(jj),:) = l(jj,:)*dx+d1-dx/2;
            end
        end
    else % if no new param vector, only keep the initial one
        X = P0.pts(:,ii);
        nepsi = P0.epsi(:,ii);
    end
    
    P.pts = [P.pts, X]; % we cant know a priori the size of P.pts
    P.epsi = [P.epsi, nepsi];
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
P.traj_ref = zeros(1,size(P.pts,2));

% check if some already computed trajectories are valid for generated param
% set
[Ppts,~,i_pts] = unique(P.pts(1:P.DimP,:)','rows','stable'); % got unique param vector for P
[~,iP,iP0] = intersect(Ppts,P0.pts(1:P.DimP,logical(P0.traj_ref))','rows','stable'); % get param vector common between computed P0 and (unique) P

for ii = 1:numel(iP) % for all param vector common between P and those computed in P0
    P.traj_ref(i_pts==iP(ii)) = ii; % affect the ith traj to all param vector equals to the one in P0
    P.traj(ii) = P0.traj(P0.traj_ref(iP0(ii))); % and fill the ith traj in P
end

% Then define which are remaining trajectories to compute
[~,P.traj_to_compute] = unique(P.pts(1:P.DimP,:)','rows','first');
P.traj_to_compute = setdiff(P.traj_to_compute,find(P.traj_ref~=0)); % don't keep those already computed
P.traj_to_compute = sort(P.traj_to_compute)'; % set it in a line shape

end
