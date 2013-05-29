function P = SobolRefine(P, nb, varargin)
% SOBOLREFINE  Sample quasi-uniformly a parameter set using Sobol sequence
%
% Credit:  John Burkardt, 2006
%
% Synopsis:  P = SobolRefine(P, nb [ , step ] [ , 'striclyInside' ] )
%
% Input:
%  - P    : the parameter set to refine. It may contains one or many
%           set of parameter values
%  - nb   : the number of new set of parameter to generate for each
%           parameter set in P. If nb is lower or equal to one, P itself is
%           answered.
%  - step : number of initial generated values by sobol generator to skip
%           (optional, default=0).
%  - striclyInside : if the keyword 'striclyInside' is provided, the new
%                    parameter set are created such that all new boxes are
%                    contained in the initial box. Otherwise, the new boxes
%                    are created such that the center of the new boxes are
%                    in the initial one, but some new boxes may be such
%                    that Pnew.pts - Pnew.epsi < Pold.pts - Pold.epsi or
%                    Pnew.pts + Pnew.epsi > Pold.pts + Pold.epsi.
%                    (optional, not set by default)
%
% Output:
%  - P : the new parameter set
%
% Example (Lorentz84):
%
%   CreateSystem;
%   P = CreateParamSet(Sys,{'F','G','a','b'},[1,100;0,5;.24,.26;2.9,3.1]);
%   
%   Pr = SobolRefine(P, 4, 'striclyInside');
%   Pr = SobolRefine(Pr, 200, 'striclyInside'); % Sample with 800 points
%   
%   Pr2 = SobolRefine(P, 4); % 4 parameters sets
%   Pr2 = SobolRefine(Pr2, 200) % also 800 parameter sets
%
%   Pr3 = SobolRefine(P, 800, 'striclyInside');
%
%   Pr4 = SobolRefine(P, 800);
%
%   figure;
%   subplot(2,3,1);
%   SplotBoxPts(P); % Parameter set before sampling
%   SplotPts(Pr);   % plots the generated points using striclyInside
%   subplot(2,3,2);
%   SplotBoxPts(P);
%   SplotBoxPts(Pr);
%   subplot(2,3,3);
%   SplotBoxPts(P);
%   SplotPts(Pr3);
%   subplot(2,3,4);
%   SplotBoxPts(P);   % plots the generated points without using striclyInside
%   SplotPts(Pr2);
%   subplot(2,3,5);
%   SplotBoxPts(P);
%   SplotBoxPts(Pr2);
%   subplot(2,3,6);
%   SplotBoxPts(P);
%   SplotPts(Pr4);
%
%See also QuasiRefine Refine RandomLogRefine
%

if(nb<=1)
    return ;
end

if(nargin==2)
    step=0;
    striclyInside = false;
elseif(nargin==3)
    if ischar(varargin{1})
        striclyInside = true;
        step = 0;
    else
        striclyInside = false;
        step = varargin{1};
    end
else
    step = varargin{1};
    striclyInside = true;
end


dim_num = numel(P.dim);

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
if(striclyInside)
    width = 2*(old_epsi - new_epsi);
    mini = P.pts(P.dim,:) - (old_epsi - new_epsi);
else
    width = 2 * old_epsi;
    mini = P.pts(P.dim,:) - old_epsi;
end
P.pts = repmat(P.pts,[1 nb]);
P.pts(P.dim,:) = repmat(width,[1 nb]).*r+repmat(mini,[1 nb]);


if isfield(P,'selected')
    P.selected = zeros(1, size(P.pts,2));
end

end
