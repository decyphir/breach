function [PRLog] = CreateRandomLogParamSets(Sys, params, ranges, N, varargin)
%CREATERANDOMLOGPARAMSETS Create a logarithmic random sampling of parameters.
% The unknown parameters are defined by the argument params. If the lower
% limit of an interval is zero, it is replaced by the absolute tolerance of
% the system. There should not be a negative range limit. It must not be
% two (or more) times the same parameter in the params argument.
%
% Synopsis: PRLog = CreateRandomLogParamSets(Sys, params, ranges, N [ , 'strictlyInside' ] )
%
% Inputs:
%  - Sys    : The considered system
%  - params : List of varying parameters
%  - ranges : Parameters ranges
%  - N      : Number of random generated points
%  - 'strictlyInside' : (Optional, default=not set) If set, all generated
%                       parameter sets are strictly within the ranges.
%                       Otherwise, only the center of the generated
%                       parameter sets are within the provided ranges.
%
% Output:
%  - PRLog : A random logarithmic sampling of N points
%
% Example (Lorentz84):
%    CreateSystem;
%    PLog = CreateRandomLogParamSets(Sys,{'a','b'},[1.0e-3,1.0e2;3,5],50);
%            % PLog is a parameter set containing 50 parameter values such
%            % that a is logarithmicly distributed in [1.0e-3 , 1.0e2]
%            % and b in [3 , 5]
%    SplotBoxPts(PLog);
%
%See also CreateParamSet RandomLogRefine SAddUncertainParam
%

% we check if all range limits are strictly positive
if ~isempty(find(sign(ranges)<0,1))
    %return;
    error('CreateRandomLogParamSets:rangeBound',...
            'Range limits must be positive or null.');
end

if ischar(params)
    params = {params};
end

nbParam = numel(params);

% check for duplicates in params
if(nbParam~=numel(unique(params)))
    error('CreateRandomLogParamSets:duplicateParam',...
            ['When refining, it must not be more than one time a '...
            'parameter']);
end

strictlyInside = nargin==5 && strcmpi(varargin{1},'strictlyinside');


% if there are null range limits, we replace them by the AbsTol
ranges(ranges==0) = Sys.CVodesOptions.AbsTol;

ranges = log10(ranges); %from normal scale to log scale

PRLog = CreateParamSet(Sys, params, ranges);
if(strictlyInside)
    PRLog = QuasiRefine(PRLog,N,'strictlyInside');
else
    PRLog = QuasiRefine(PRLog,N);
end

for ii=1:nbParam % we rescale each param, in the order of epsi
    epsi = PRLog.epsi(ii,:);
    value = PRLog.pts(PRLog.dim(ii),:);
    
    vinf = value-epsi;
    vsup = value+epsi;

    vinf = 10.^vinf; %lowest new value on normal scale
    vsup = 10.^vsup; %highest new value on normal scale
    epsi = (vsup-vinf)/2; %new epsi on normal scale
    value = vinf + epsi; %compute new value on normal scale    
    
    PRLog.epsi(ii,:) = epsi;
    PRLog.pts(PRLog.dim(ii),:) = value;
end

end
