function Sr = HaltLogRefine(S, N, l,u)
%HaltLogRefine  sample a parameter set uniformly on a logarithmic scale
%
%  Ex:
%      CreateSystem;
%      S  = CreateParamSet(Sys, [1 2], [0 2;0 2]');
%      Sr = HaltLogRefine(S, 100, -2,2); % log sampling 100 samples between .01 and 100
%

n = numel(S.dim);

if(size(P.pts,2)~=1)
    error('Not implemented for multiple points. P should have exactly one point.')
end

if(numel(nargin)==1)
    l=-2; u=2; % by default, goes between 1/00 and 100 times the nominal value
end

pts_val = (u+l)/2;
epsi_val = (u-l)/2;

Stmp = SCreate(pts_val*ones(n,1),epsi_val*ones(n,1));
Stmp = HaltonRefine(Stmp,N);

Stmp.pts = 10.^(Stmp.pts).*repmat(S.pts(S.dim,:), [1 N]);

Sr = S;
Sr.pts = repmat(S.pts,[1 N]);

Sr.pts(S.dim, :) = Stmp.pts;
Sr.epsi = repmat(S.epsi./N,[1 n]);

end
