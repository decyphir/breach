function P = SobolRefine(P, nb, varargin)
%SOBOLREFINE samples quasi-uniformly a parameter set using Sobol sequence.
% This function should not be used directly by through QuasiRefine function
% as it does not manage the traj_ref and traj_to_compute fields.
%
% Credit:  John Burkardt, 2006
%
% Synopsis:  P = SobolRefine(P, nb[ , step][, 'striclyInside'])
%
% Inputs:
%  - P    : the parameter set to refine. It may contain many parameter
%           vectors
%  - nb   : the number of new parameter vectors generated for each
%           parameter vector of P. If nb is lower or equal to one, P is not
%           modified.
%  - step : (Optional, default=0) number of initial generated values by
%           sobol generator to skip
%  - striclyInside : if the keyword 'striclyInside' is provided, the new
%                    parameter vectors are created such that all new boxes
%                    are contained in the initial box. Otherwise, the new
%                    boxes are created such that the center of the new
%                    boxes are in the initial one, but some new boxes may
%                    be such that Pnew.pts - Pnew.epsi < Pold.pts -
%                    Pold.epsi or Pnew.pts + Pnew.epsi > Pold.pts +
%                    Pold.epsi. (optional, not set by default)
%
% Output:
%  - P : the new parameter set
%
% Example (Lorentz84):
%   CreateSystem
%   P = CreateParamSet(Sys,{'F','G','a','b'},[1,100;0,5;.24,.26;2.9,3.1]);
%   
%   Pr = QuasiRefine(P, 4, 'striclyInside', 'Sobol');
%   Pr = QuasiRefine(Pr, 100, 'strictlyInside', 'Sobol'); % Sample with 400 points
%   
%   Pr2 = QuasiRefine(P, 4, 'Sobol'); % 4 parameters sets
%   Pr2 = QuasiRefine(Pr2, 100, 'Sobol'); % also 400 parameter sets
%
%   Pr3 = QuasiRefine(P, 400, 'strictlyInside', 'Sobol');
%
%   Pr4 = QuasiRefine(P, 400, 'Sobol');
%
%   figure
%   subplot(2,3,1);
%   SplotBoxPts(P,[],[],[],'g'); % Parameter set before sampling
%   SplotPts(Pr);   % plots the generated points using striclyInside
%   title('Pr : With strictlyInside');
%   subplot(2,3,2);
%   SplotBoxPts(P,[],[],[],'g');
%   SplotBoxPts(Pr);
%   title('Pr : With strictlyInside');
%   subplot(2,3,3);
%   SplotBoxPts(P,[],[],[],'g');
%   SplotPts(Pr3);
%   title('Pr3 : With strictlyInside');
%   
%   subplot(2,3,4);
%   SplotBoxPts(P,[],[],[],'g');   % plots the generated points without using striclyInside
%   SplotPts(Pr2);
%   title('Pr2 : Without strictlyInside');
%   subplot(2,3,5);
%   SplotBoxPts(P,[],[],[],'g');
%   SplotBoxPts(Pr2);
%   title('Pr2 : Without strictlyInside');
%   subplot(2,3,6);
%   SplotBoxPts(P,[],[],[],'g');
%   SplotPts(Pr4);
%   title('Pr4 : Without strictlyInside');
%   
%   pos = get(gcf, 'Position');
%   set(gcf, 'Position', [pos(1), pos(2), pos(3)*1.5, pos(4)*1.5]);
% 
%See also QuasiRefine Refine RandomLogRefine LogNRefine
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
%P.epsi = kron(P.epsi, ones(1,size(P.pts,2)))/nb;

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

P = Preset_traj_ref(P);

end
