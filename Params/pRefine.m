function P = pRefine(P, p, r)
%
% PREFINE refines the parameter set P by picking r initial points for the
% Morris global sensitivity measure. These points are randomly chosen in
% the grid with p levels (See Saltelli's Books, chapter 3 or 4). We use the
% suggested value for Delta: p/(2(p-1))
%
% Usage:  Pr = pRefine(P, p, r)
%
% Input:
%   - P   the initial parameter set
%   - p   indicates in how many the interval of each parameter is splitted.
%         To get the sensitivity, we compute between a simulation and the
%         simulation obtained by adding or removing p/(2(p+1)) to one
%         parameter
%   - r   the number of initial points for paths to generate
%
% Output:
%   - Pr  the refined parameter set

% The method consists to generate r paths of n+1 (or n?) parameter sets
% (where n is the number of parameters), such that :
%  1/ only one parameter differs between two consecutive steps ;
%  2/ this parameter differs by p/(2(p+1)) times the size of the interval
% of the changing parameter
%  3/ all the parameters change once along the path

% define the admissible grid and pick points

n = numel(P.dim);           % dimension
delta = p/(2*(p-1));
ngrid = floor(p*(1-delta));

% we define the space for the initial points of the paths, in the unity box
Stmp.dim = 1:n;
Stmp.pts = (ngrid+1)/2*ones(n,1);
Stmp.epsi = (ngrid+1)/2*ones(n,1);
Stmp.DimP = n;
Stmp.DimX = n;

Stmp = QuasiRefine(Stmp,r);
Stmp.pts = floor(Stmp.pts)/(p-1);

% can be improved
S2 = Stmp;
S2.pts = [];
S2.epsi = [];
S2.D = [];

% generate the elementary effects "trajectories"
for k=1:r
    [X, D] = EE_traj(Stmp.pts(:,k), p, n);
    S2.pts =  [S2.pts X];
    S2.epsi =  [S2.epsi repmat(Stmp.epsi(:,k), [1 n+1] )];
    S2.D = [S2.D D];
    
end

% Normalize to P ranges
P.pts = repmat(P.pts, [1 size(S2.pts,2)]);
P.epsi = repmat(P.epsi, [1 size(S2.pts,2)]);
P.pts(P.dim,:) = P.pts(P.dim,:) + (2*S2.pts-1).*P.epsi;
P.D = S2.D;

end
