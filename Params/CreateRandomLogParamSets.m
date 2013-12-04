function [PRLog] = CreateRandomLogParamSets(Sys, params, ranges, N, varargin)
%CREATERANDOMLOGPARAMSETS creates a logarithmic random sampling of
% parameters. The unknown parameters are defined by the argument params. If
% the lower limit of an interval is zero, it is replaced by the absolute
% tolerance of the system. There should not be a negative range limit. It
% cannot be two (or more) times the same parameter in the params argument.
% 
% Synopsis: PRLog = CreateRandomLogParamSets(Sys, params, ranges, N[, 'strictlyInside'])
% 
% Inputs:
%  - Sys    : The considered system
%  - params : List of names or indexes of varying parameters
%  - ranges : array of size numel(params) x 2 describing the parameters
%             ranges
%  - N      : Number of random generated points
%  - 'strictlyInside' : (Optional, default=not set) If set, all generated
%                       parameter sets are strictly within the ranges.
%                       Otherwise, only the center of the generated
%                       parameter sets are within the provided ranges.
% 
% Output:
%  - PRLog : A random logarithmic sampling of N parameter vectors
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
if numel(Sys.CVodesOptions.AbsTol)==1
    ranges(ranges==0) = Sys.CVodesOptions.AbsTol;
else
    ranges(ranges==0) = Sys.CVodesOptions.AbsTol(ranges==0);
end

ranges = log10(ranges); %from normal scale to log scale

PRLog = CreateParamSet(Sys, params, ranges);
if(strictlyInside)
    PRLog = QuasiRefine(PRLog,N,'strictlyInside');
else
    PRLog = QuasiRefine(PRLog,N);
end

% Translate back to linear scale, using the midpoint and the upper bound to
% compute pts and epsi (using lower bound and upper bound may lead to pts
% beyond intervals)
val = PRLog.pts(PRLog.dim,:);
PRLog.pts(PRLog.dim,:) = 10.^val;
PRLog.epsi = 10.^(val+PRLog.epsi) - PRLog.pts(PRLog.dim,:); % epsi = sup - value

end
