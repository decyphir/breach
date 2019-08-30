function P = Refine(P0, delta, sample_bndary)
%REFINE generates a grid points in a N-dimensional parameter set.
% 
% Synopsis: P = Refine(P0, delta)
% 
% Inputs:
%   - P0    : is a parameter set which may contain many parameter vectors
%   - delta : is a vector of dimension numel(P0.dim) x 1. It also may be a
%             scalar. In this case, it is interpreted as a column vector
%             with all components equal to its value. The ith line
%             corresponds to the ith parameter in P0.dim. Its values can be
%             either integers or reals. If they are integers, P has
%             delta(i) points in dimension i. If any is real (ie:
%             any(floor(delta)~=delta) is true), then delta is interpreted
%             as a distance and P is divided into points that are at
%             distance delta from one another. If any value in delta in
%             lower or equal to 0, P0 is answered.
% 
% Output:
%  - P : is the generated parameter set
% 
% Example (lorentz84):
%   %First, create some parameter set, 2 dimensions ranging from 0 to 2:
%   
%   CreateSystem;
%   P0 = CreateParamSet(Sys, [1, 2], [0, 2; 0, 2]);
%   
%   %Then :
%   
%   P = Refine(P0, 2); % creates a 2x2 grid in [0 2 0 2]
%   P.pts              % show the generated parameter vectors
%   
%   P = ComputeTraj(Sys,P0,1:0.1:10);
%   P = Refine(P, [3;3]);   % creates a 3x3 grid
%   P.traj_ref              % should be [0 0 0 0 1 0 0 0 0]
%   P.traj_to_compute       % should be [1 2 3 4 6 7 8 9]
%   P.traj
%   
%   P = Refine(P0, .1);    % creates a grid with resolution (.1,.1)
%   P.pts(:,[1,2,3,19,20,21,22,400])  % show some parameter vectors: the
%                                     % difference between to parameter
%                                     % vectors is 0.1
%   
%   Refine(P0, [.1;.2]) % creates a grid with resolution (.1,.2) (works
%                       % well as .1 and .2 are divisor of 2*P.epsi)
%   
%   P = Refine(P0, [1.3;1.3]); % creates a grid with "resolution (1.3,1.3)"
%   P.pts    % the difference between two consecutive values is not 1.3 as
%            % 1.3 is not a divisor of 2*P.epsi
%   
%   P = Refine(P0, [2.1;2.2]);   % returns P0:
%   all(size(P.pts)==size(P0.pts)) & all(P.pts==P0.pts)  % true
% 
%See also RandomLogRefine QuasiRefine CountRefine CreateParamSet
%SAddUncertainParam SetEpsi
%

if ~exist('sample_bndary','var')
    sample_bndary= 0;
end

% check delta
if(all(delta==1) || any(delta<=0))
    P = P0;
    return;
end

if isscalar(delta)
    delta = delta*ones(numel(P0.dim),1);
elseif(size(delta,1)==1)
%    warning('Refine:BadDeltaShape','delta must have the format numel(P0.dim) x 1. Transposed.');
    delta = delta'; % try to transpose if needed
end

if(numel(delta)>numel(P0.dim))
    warning('Refine:TooLargeDelta','The number of element in delta is higher than P0.dim. Only first elements kept.');
    delta = delta(1:numel(P0.dim));
elseif(numel(delta)<numel(P0.dim))
    error('Refine:SizeOfDelta','The number of element in delta is different than P0.dim.');
end

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
            if(deltai(jj)>1)
                d1 = xlim(jj,1);
                d2 = xlim(jj,2);
                if sample_bndary
                    dx = (d2-d1)./(deltai(jj)-1);
                    X(P0.dim(jj),:) = (l(jj,:)-1)*dx+d1;           
                else
                    dx = (d2-d1)./(deltai(jj));
                    X(P0.dim(jj),:) = l(jj,:)*dx+d1-dx/2;
                end
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

if isfield(P0,'time_mult')
    P.time_mult = P0.time_mult;
end

P.DimX = P0.DimX;
P.DimP = P0.DimP;

P = Preset_traj_ref(P);


end
