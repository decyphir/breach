function Y = GetTrajValues(Pf, iX,  t)
% GETTRAJVALUES  extract values of a variable evolution in a set of
% trajectories
%
% Synopsis: Y = GetTrajValues(Pf,iX,t)
%
% Input:
%  -  Pf : parameter set with pre-computed trajectories
%  -  iX : index or name of the required variable
%  -  t  : time instant(s) when to extract values. Can be a scalar or an
%          array of time instants
%
% Output:
%  - Y : the value of iX at time instant(s) t. Each column contains the
%        value of iX for the different trajectories at a timestep
%
% Example (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys,{'a','b'},[0.2,0.3;3,5]);
%   P = Refine(P,2);
%   P = ComputeTraj(Sys,P,0:0.1:10);
%   Y = GetTrajValues(P,'x0',2:1:8);
%

if (~isfield(Pf,'traj'))
    error('GetTrajValues:NoTrajField','Compute trajectories first')
end

if ischar(iX) || iscell(iX)
    iX = FindParam(Pf, iX);
end

%iX = iX(iX<=Pf.DimX); % NM: we keep only the variables, not mandatory if
%GetTrajValues is correctly used

if (~exist('t','var'))
    t = Pf.traj{1}.time;
end

X = cat(1, Pf.traj{1}.X); % concatenate all trajectories

X = X(iX:Pf.DimX:end,:); % keep only the evolution of iX over time

Y = interp1(Pf.traj{1}.time, X',t)';

end
