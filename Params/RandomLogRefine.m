function [P] = RandomLogRefine(P, N, minValue)
%RANDOMLOGREFINE Create a logarithmic random sampling of parameters. If P
%   contains many points, each of them is divided into N new points. All
%   ranges (ie : [value-epsi ; value+epsi] must be strictly positive. When
%   there is more than 16 ordre of magnitude (with Matlab R2012), it is
%   recommended to use CreateRandomLogParamSets or you risk to face error
%   due to null interval limit.
%
% Syntax: PRLog = RandomLogSampling(P, N [, minValue] )
%
% Inputs:
%
%    -  P         The parameter set to refine. If the number of uncertain
%                 parameters higher than 40, the function may thrown an
%                 error.
%    -  N         Number of random generated points. If N is lower or equal
%                 to one, nothing is done.
%    -  minValue  If there are parameters with intervals such that
%                 "pts-epsi<=0", if minValue exists, the lower possible
%                 value for this parameter is set to minValue, otherwize,
%                 an error is thrown
%
% Outputs:
%
%    -  PRLog    A random logarithmic sampling of N points
%
% Example (lorentz84):
%
%   CreateSystem;
%   P = CreateParamSet(Sys);
%   PRLog = RandomLogRefine(P,10);
%   SplotBoxPts(PRLog);    % Parameter set after sampling
%
%TODO : THIS FUNCTION DOES NOT MANAGE THE TRAJ_REF AND TRAJ_TO_COMPUTE
% FIELDS
%
%See also CreateRandomLogParamSets Refine
%

% 0/ check for input
if(N<=1)
    return;
end


% 1/ convert from unity scale to log scale
epsi = P.epsi;
val = P.pts(P.dim,:);
mini = val-epsi; %we assume inf<=sup

if any(mini<=0)
    %return;
    if exist('minValue','var')
        mini(mini<=0)=minValue;
    else
        error('RandomLogRefine:rangeBound','Range limits must be strictly positive.');
    end
end

mini = log10(mini);
maxi = log10(val+epsi);
P.pts(P.dim,:) = (maxi+mini)/2;
P.epsi = (maxi-mini)/2;


% 2/ compute quasi-refine on log scale
ii=0; % look for highest power of 2 lower than N
while(2^ii<N) % N is > 1
    ii = ii + 1;
end
P = QuasiRefine(P, N, 2^(ii-1));


% 3/ convert back from log scale to unity scale
epsi = P.epsi;
val = P.pts(P.dim,:);
mini = 10.^(val-epsi);
maxi = 10.^(val+epsi);
P.pts(P.dim,:) = (maxi+mini)/2;
P.epsi = (maxi-mini)/2;


end
