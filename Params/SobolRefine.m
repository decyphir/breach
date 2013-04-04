function P = SobolRefine(P, nb, step)
% SOBOLREFINE  Sample quasi-uniformly a parameter set using Sobol sequence
%
% Credit:  John Burkardt, 2006
%
% Synopsis:  Ph = SobolRefine(P, nb, step)
%
% Example (Lorentz84):
%
%   CreateSystem;
%   P = CreateParamSet(Sys,{'F','G'},[1,100;0,5]);
%   Ph = SobolRefine(P, 1000); % Sample with 1000 points
%   
%   Ph2 = Refine(P, 2); % 4 parameters sets
%   Ph2 = SobolRefine(Ph2, 250) % also 1000 parameter sets
%
%   SplotBoxPts(P); % Parameter set before sampling
%   SplotPts(Ph2);   % plots the generated points
%
%See also QuasiRefine Refine RandomLogRefine
%

if(nb<=1)
    return ;
end

dim_num = numel(P.dim);

if (nargin==2)
    step=0;
end

%generate random values between 0 and 1
r = i4_sobol_generate(dim_num, nb, step);
r = kron(r, ones(1,size(P.pts,2)));

width = 2*P.epsi;   % NM : TODO make modifications to avoid parameter sets out of the initial interval
inf = P.pts(P.dim,:)-P.epsi;   % NM : TODO make modifications to avoid parameter sets out of the initial interval
P.pts = repmat(P.pts,[1 nb]);
P.pts(P.dim,:) = repmat(width,[1 nb]).*r+repmat(inf,[1 nb]);

P.epsi = repmat(P.epsi,[1 nb])/(floor(nb^(1/dim_num)));

% Wrong : the dimension of epsi differs from the dimension of pts
%Ph.epsi = kron(Ph.epsi,ones(1,size(Ph.pts,2)))/(floor(nb^(1/dim_num)));

% used to avoid superposition of square when selectionned and shown in 2
% dimension ; it is not really correct because space is "lost"
%Sh.epsi = kron(Sh.epsi, ones(1,size(Sh.pts,2)))/nb;

if (isfield(P,'selected'))
    P.selected = zeros(1, size(P.pts,2));
end

end
