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

% A AMELIORER :
%
% Utiliser la génération de nombre aléatoire type quasi monte-carlo au lieu
% du rand proposé par matlab.

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

PRLog = CreateSampling(Sys, params); %pas besoin des intervalles
PRLog = Refine(PRLog, [N,ones(1,nbParam-1)]);  %On découpe en N points

% if there are null range limits, we replace them by the AbsTol
posNull = ranges==0;
ranges(posNull) = Sys.CVodesOptions.AbsTol;

ranges = log10(ranges);

for i=1:nbParam
    %Pour chaque param, on défini une position aléatoire
    epsi = (ranges(i,2)-ranges(i,1))/2; %epsi de l'intervalle initial en log
    r = ranges(i,1) + 2*epsi*rand(1,N); %on tire au hasard sur l'échelle log
    epsi = epsi/(N^(1/nbParam)); %epsi des nouveaux intervalles en log

    % possibility 1 : define epsi, using the lower bound
    inf = r - epsi; %valeur minimal des nouveaux intervalles en log
    value = 10.^r; %valeur des nouveaux points en normal
    inf = 10.^inf; %valeur minimal des nouveaux intervalles en normal
    epsi = value - inf; %epsi des nouveaux intervalles en normal

    %possiblity 2 : define epsi using the higher bound
%    sup = r + epsi; %valeur maximale des nouveaux intervalles en log
%    value = 10.^r; %valeur des nouveaux points en normal
%    sup = 10.^sup; %valeur maximale des nouveaux intervalles en normal
%    epsi = sup - value; %epsi des nouveaux intervalles en normal

    %possibility 3 : define value using the lower and higher bounds
%    inf = r - epsi; %valeur minimal des nouveaux intervalles en log
%    sup = r + epsi; %valeur maximale des nouveaux intervalles en log
%    inf = 10.^inf; %valeur minimal des nouveaux intervalles en normal
%    sup = 10.^sup; %valeur maximale des nouveaux intervalles en normal
%    epsi = (sup-inf)/2; %epsi des nouveaux intervalles en normal
%    value = inf + epsi; %valeur des nouveaux points en normal
    
    PRLog = SetParam(PRLog,params(i),value);
    position = FindParam(Sys,params{i}); %on cherche la position des epsi
    position = PRLog.dim==position;
    PRLog.epsi(position,:) = epsi;
end

end
