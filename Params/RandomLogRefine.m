function P = RandomLogRefine(P, N, minValue)
%RANDOMLOGREFINE creates a logarithmic random sampling of parameters. If P
% contains many parameter vectors, each of them is divided into N new
% parameter vectors. All ranges (ie : [value-epsi ; value+epsi] must be
% strictly positive. When there is more than 16 ordres of magnitude (with
% Matlab R2012) between the lower bound (value-epsi) and the upper bound
% (value+epsi), it is recommended to use CreateRandomLogParamSets or you
% risk to face error due to interval limit equals to zero.
% 
% Synopsis: P = RandomLogSampling(P, N[, minValue])
% 
% Inputs:
%  - P        : The parameter set to refine. It may contain many parameter
%               vectors
%  - N        : Number of random generated points for each initial
%               parameter vector. If N is lower or equal to one, nothing
%               is done. If N is not an integer, it is rounded toward 0.
%  - minValue : (Optional) if there are parameters with intervals such
%               that "pts-epsi<=0", if minValue exists, the lower possible
%               value for this parameter is set to minValue, otherwise, an
%               error is thrown.
% 
% Outputs:
%  - P : A random logarithmic sampling of N parameter vectors time the
%        initial number of parameter vectors. The central value of each
%        new parameter vector is included in the initial set of
%        hyper-rectangle described by the parameter vectors. Nevertheless,
%        the new hyper-rectangle may extand beyond the initial
%        hyper-rectangles.
% 
% Example (lorentz84):
%   CreateSystem
%   P = CreateParamSet(Sys,{'x0','x1'},[1,5;2,10]);
%   PRLog = RandomLogRefine(P,10);
%   SplotBoxPts(PRLog);    % Parameter set after sampling
% 
%See also CreateRandomLogParamSets Refine
%

% 0/ check for input
if(N<=1)
    return;
end

N = floor(N);


% 1/ convert from linear scale to log scale
epsi = P.epsi;
val = P.pts(P.dim,:);
mini = val-epsi;

if any(mini<=0)
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


% 3/ convert back from log scale to linear scale
% we cannot compute the lower and higher bound of the new interval and
% deduce pts and epsi from them because pts = (10^lower+10^higher)/2 may be
% higher than the initial upper bound
val = P.pts(P.dim,:);
P.pts(P.dim,:) = 10.^val;
P.epsi = 10.^(val+P.epsi) - P.pts(P.dim,:); % epsi = sup - value

end
