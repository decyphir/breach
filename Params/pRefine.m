function P = pRefine(P0, p, r, step)
% PREFINE refines the parameter set P by picking r initial points for the
% Morris global sensitivity measure. These points are randomly chosen in
% the grid with p levels (See Saltelli's Books, chapter 3 or 4). We use the
% suggested value for Delta: p/(2(p-1))
%
% Usage:  Pr = pRefine(P, p,r,step)
%
% Input:
%   - P   the initial parameter set. Only the first set of value of P will
%         be considered. The parameter not in P.dim are equals in all the
%         generated parameter sets. There value is provided by the first
%         set of value in P.
%   - p   indicates in how many the interval of each parameter is splitted.
%         To get the sensitivity, we compute between a simulation and the
%         simulation obtained by adding or removing p/(2(p+1)) to one
%         parameter
%   - r   the number of initial points for paths to generate
%   - step seed for quasi random sampling
% Output:
%   - Pr  the refined parameter set containing (k+1)*r parameter sets,
%         where k=numel(P.dim). A field D is set in Pr, which indicates
%         which parameter changes of value between two consecutive
%         trajectories. The size of D is k x k*r.
%
% See also SPropSensi, EE_traj, EEffects
%

% The method consists to generate r paths of k+1 parameter sets
% (where k is the number of parameters), such that :
%  1/ only one parameter differs between two consecutive steps ;
%  2/ this parameter differs by p/(2(p+1)) times the size of the interval
% of the changing parameter
%  3/ all the parameters change exactly once along the path

% define the admissible grid and pick points

if ~exist('step', 'var')
    step = 1;
end

rng(step, 'twister');
P=SPurge(P0);

k = numel(P0.dim);           % dimension

delta = p/(2*(p-1));
ngrid = floor(p*(1-delta));

% we define the space for the initial points of the paths, in the unity box
Stmp.dim = 1:k;
Stmp.pts = (ngrid+1)/2*ones(k,1);
Stmp.epsi = (ngrid+1)/2*ones(k,1);
Stmp.DimP = k;
Stmp.DimX = k;

Stmp = QuasiRefine(Stmp,r, step);
Stmp.pts = floor(Stmp.pts)/(p-1);

%intialization to avoid growth when in loop
S2 = Stmp;
S2.pts = zeros(k,(k+1)*r);
S2.epsi = zeros(k,(k+1)*r);
S2.D = zeros(k,k*r);

% generate the elementary effects "trajectories"
for ii=1:r
    [X, D] = EE_traj(Stmp.pts(:,ii), p, k);
    S2.pts(:,(ii-1)*(k+1)+1:ii*(k+1)) = X;  % define S2.pts with r bloc of k x (k+1)
    S2.epsi(:,(ii-1)*(k+1)+1:ii*(k+1)) = repmat(Stmp.epsi(:,ii), [1,k+1]);
    S2.D(:,(ii-1)*k+1:ii*k) = D;  % fill S2.D with r bloc of k x k
    
end

% Normalize to P ranges
P.dim = P0.dim;
P.pts = repmat(P0.pts(:,1), [1 size(S2.pts,2)]);
P.epsi = repmat(P0.epsi(:,1), [1 size(S2.pts,2)]);
P.pts(P.dim,:) = P.pts(P.dim,:) + (2*S2.pts-1).*P.epsi;
P.epsi = P.epsi*delta;   % <-- NOT SURE OF THAT
P.D = S2.D;

if isfield(P0, 'traj')
    P = Pimport_traj(P, P0);
else
    P = Preset_traj_ref(P);
end

end
