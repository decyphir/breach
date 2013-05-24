function P = SobolRefine(P, nb, step)
% SOBOLREFINE  Sample quasi-uniformly a parameter set using Sobol sequence
%
% Credit:  John Burkardt, 2006
%
% Synopsis:  P = SobolRefine(P, nb, step)
%
% Input:
%  - P    : the parameter set to refine. It may contains one or many
%           set of parameter values
%  - nb   : the number of new set of parameter to generate for each
%           parameter set in P.
%  - step : number of initial generated values by sobol generator to skip.
%
% Output:
%  - P : 
%
% Example (Lorentz84):
%
%   CreateSystem;
%   P = CreateParamSet(Sys,{'F','G'},[1,100;0,5]);
%   Pr = SobolRefine(P, 1000); % Sample with 1000 points
%   
%   Pr2 = Refine(P, 2); % 4 parameters sets
%   Pr2 = SobolRefine(Pr2, 250) % also 1000 parameter sets
%
%   SplotBoxPts(P); % Parameter set before sampling
%   SplotPts(Pr2);   % plots the generated points
%
%See also QuasiRefine Refine RandomLogRefine
%

if(nb<=1)
    return ;
end

dim_num = numel(P.dim);

if(nargin==2)
    step=0;
end

%generate random values between 0 and 1
r = i4_sobol_generate(dim_num, nb, step); % generate nb point of dim_num dimensions
r = kron(r, ones(1,size(P.pts,2)));

old_epsi = P.epsi;
new_epsi = P.epsi/(nb^(1/dim_num));
P.epsi = repmat(new_epsi,[1 nb]);
% old version
%P.epsi = repmat(P.epsi,[1 nb])/(floor(nb^(1/dim_num)));
% used to avoid superposition of square when selectionned and shown in 2
% dimension ; it is not really correct because space is "lost"
%Sh.epsi = kron(Sh.epsi, ones(1,size(Sh.pts,2)))/nb;

%we scale the random generated value to make then fit the intervals. New
% generated parameter sets are strictly inside the initial one.
width = 2*(old_epsi - new_epsi);
inf = P.pts(P.dim,:) - (old_epsi - new_epsi);
P.pts = repmat(P.pts,[1 nb]);
P.pts(P.dim,:) = repmat(width,[1 nb]).*r+repmat(inf,[1 nb]);


if isfield(P,'selected')
    P.selected = zeros(1, size(P.pts,2));
end

end
