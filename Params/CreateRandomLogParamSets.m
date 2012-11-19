function [PRLog] = CreateRandomLogParamSets(Sys, params, ranges, N)
%CREATERANDOMLOGPARAMSETS Create a logarithmic random sampling of parameters.
%   The unknown parameters are defined by the argument params. If the lower
%   limit of an interval is zero, it is replaced by the absolute tolerance
%   of the system. There should not be a negative range limit.
%   It must not be two (or more) times the same parameter in the params
%   argument.
%
%   Syntax: PRLog = RandomLogSampling(Sys, params, intervals, N)
%
%   Inputs:
%
%    -  Sys       The considered system
%    -  params    List of varying parameters
%    -  intervals Parameters ranges
%    -  N         Number of random generated points
%
%   Outputs:
%
%    -  PRLog    A random logarithmic sampling of N points
%
%	Example:
%
%   RandomLogSampling(Sys,{'k1','k2'},[5.0e-9,2.0e-4;2,3],50) returns a
%   parameter set containing 50 parameter values such that k1 is
%   logarithmicly distributed in [5.0e-9 , 2.0e-4] and 'k2' in [2 , 3]
%

% TO IMPROVE :
%
% Use quasi monte-carlo generation to generate random number instead of
% using the rand operator provided by Matlab.

% we check if all range limits are strictly positive
if(~isempty(find(sign(ranges)<0,1)))
    %return;
    error('RandomLogSampling:rangeBound','Range limits must be positive or null.');
end

if(ischar(params))
    params={params};
end

nbParam = numel(params);

% check for duplicates in params
if(nbParam~=numel(unique(params)))
    %return ;
    error('RandomLogSampling:duplicateParam',['When refining, it must '...
            'not be more than one time a parameter']);
end

PRLog = CreateSampling(Sys, params); %dont need to provide intervalles
PRLog = Refine(PRLog, [N,ones(1,nbParam-1)]);  %We split into N points

% if there are null range limits, we replace them by the AbsTol
ranges(ranges==0) = Sys.CVodesOptions.AbsTol;

ranges = log10(ranges);

for i=1:nbParam
    %For each parameter, we define a random position
    epsi = (ranges(i,2)-ranges(i,1))/2; %initial epsi on log scale
    r = ranges(i,1) + 2*epsi*rand(1,N); %random value on log scale
    epsi = epsi/(N^(1/nbParam)); %new epsi on log scale

    % possibility 1 : define epsi, using the lower bound
    inf = r - epsi; %lowest value of new intervalles on log scale
    value = 10.^r; %value of new points on normal scale
    inf = 10.^inf; %minimal value of the intervalles on normal scale
    epsi = value - inf; %compute new epsi in normal scale

    %possiblity 2 : define epsi using the higher bound
%    sup = r + epsi; %highest value for new intervalles on log scale
%    value = 10.^r; %new points on normal scale
%    sup = 10.^sup; %highest value on normal scale
%    epsi = sup - value; %compute new epsi on normal scale

    %possibility 3 : define value using the lower and higher bounds
%    inf = r - epsi; %lowest new value on log scale
%    sup = r + epsi; %highest new value on log scale
%    inf = 10.^inf; %lowest new value on normal scale
%    sup = 10.^sup; %highest new value on normal scale
%    epsi = (sup-inf)/2; %new epsi on normal scale
%    value = inf + epsi; %compute new value on normal scale

    PRLog = SetParam(PRLog,params(i),value);
    position = FindParam(Sys,params{i}); %look for epsi index in ParamList
    position = PRLog.dim==position; %then for the line in the epsi array
    PRLog.epsi(position,:) = epsi;
end

end
