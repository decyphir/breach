function Ph = SobolRefine(P, nb, step)
% SOBOLREFINE  Sample quasi-uniformly a parameter set using Sobol sequence
%
% Synopsis:  Ph = SobolRefine(P, nb)
%
% Example:
%
%   CreateSystem;
%   P = CreateSampling(Sys); % Create default parameter set for system Sys
%   Ph = SobolRefine(P, 1000); % Sample with 1000 points
%
%   SplotBoxPts(P); % Parameter set before sampling
%   SplotPts(Ph);   % plots the generated points
%
% Credit:  John Burkardt, 2006
%

dim_num = numel(P.dim);

if (nargin==2)
    step=0;
end

%generate random values between 0 and 1
r = i4_sobol_generate(dim_num, nb, step);
r = kron(r, ones(1,size(P.pts,2)));

A = 2*P.epsi;
a = P.pts(P.dim,:)-P.epsi;
Ph = P;
Ph.pts = repmat(P.pts,[1 nb]);
Ph.pts(P.dim,:) = repmat(A,[1 nb]).*r+repmat(a,[1 nb]);

% Why the following line?
%Sh.epsi = repmat(Sh.epsi,[1 nb])/(floor(nb^(1/dim_num)));

Ph.epsi = kron(Ph.epsi,ones(1,size(Ph.pts,2)))/(floor(nb^(1/dim_num)));

% used to avoid superposition of square when selectionned and shown in 2
% dimension ; it is not really correct because space is "lost"
%Sh.epsi = kron(Sh.epsi, ones(1,size(Sh.pts,2)))/nb;

if (isfield(P,'selected'))
    Ph.selected = zeros(1, size(P.pts,2));
end

end
