function [PRLog] = RandomLogRefine(P, N, minValue)
%RANDOMLOGREFINE Create a logarithmic random sampling of parameters. If P
%   contains many points, each of them is divided into N new points. All
%   ranges (ie : [value-epsi ; value+epsi] must be strictly positive. When
%   there is more than 16 ordre of magnitude (with Matlab R2012), it is
%   recommended to use CreateRandomLogParamSets or you risk to face error
%   due to null interval limit.
%
%   Syntax: PRLog = RandomLogSampling(P, N [,minValue] )
%
%   Inputs:
%
%    -  P         The parameter set to refine
%    -  N         Number of random generated points
%    -  minValue  If there are parameters with intervals such that
%                 "pts-epsi<=0", if minValue exists, the lower possible
%                 value for this parameter is set to minValue, otherwize, a
%                 error is thrown
%
%   Outputs:
%
%    -  PRLog    A random logarithmic sampling of N points
%
%	Example (lorentz84):
%
%   CreateSystem;
%   P = CreateParamSet(Sys);
%   PRLog = RandomLogRefine(P,10);
%   SplotBoxPts(PRLog);    % Parameter set after sampling
%
% SEE ALSO CREATERANDOMLOGPARAMSETS
%


PRLog = P;


% convert from unity scale to log scale
for i=1:numel(P.dim)
    epsi = P.epsi(i,:);
    val = P.pts(P.dim(i),:);
    inf = val-epsi; %we assume inf<=sup
    
    if any(inf<=0)
        %return;
        if exist('minValue','var')
            inf(inf<=0)=minValue;
        else
            error('RandomLogSampling:rangeBound','Range limits must be strictly positive.');
        end
    end
    
    inf = log10(inf);
    sup = log10(val+epsi);
    PRLog.pts(PRLog.dim(i),:) = (sup+inf)/2;
    PRLog.epsi(i,:) = (sup-inf)/2;
end


% compute quasi-refine on log scale
PRLog = QuasiRefine(PRLog, N);


% convert back from log scale to unity scale
for i=1:numel(PRLog.dim)
    epsi = PRLog.epsi(i,:);
    val = PRLog.pts(PRLog.dim(i),:);
    inf = 10.^(val-epsi);
    sup = 10.^(val+epsi);
    PRLog.pts(PRLog.dim(i),:) = (sup+inf)/2;
    PRLog.epsi(i,:) = (sup-inf)/2;
end


end
