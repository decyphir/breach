function [PRLog] = CreateRandomLogParamSets(Sys, params, ranges, N)
%CREATERANDOMLOGPARAMSETS Create a logarithmic random sampling of parameters.
%   The unknown parameters are defined by the argument params. If the lower
%   limit of an interval is zero, it is replaced by the absolute tolerance
%   of the system. There should not be a negative range limit.
%   It must not be two (or more) times the same parameter in the params
%   argument.
%
% Syntax: PRLog = CreateRandomLogParamSets(Sys, params, ranges, N)
%
% Inputs:
%
%    -  Sys       The considered system
%    -  params    List of varying parameters
%    -  ranges    Parameters ranges
%    -  N         Number of random generated points
%
% Outputs:
%
%    -  PRLog    A random logarithmic sampling of N points
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
    params={params};
end

nbParam = numel(params);

% check for duplicates in params
if(nbParam~=numel(unique(params)))
    %return ;
    error('CreateRandomLogParamSets:duplicateParam',...
            ['When refining, it must not be more than one time a '...
            'parameter']);
end

% if there are null range limits, we replace them by the AbsTol
ranges(ranges==0) = Sys.CVodesOptions.AbsTol;

ranges = log10(ranges); %from normal scale to log scale

PRLog = CreateParamSet(Sys, params, ranges);
PRLog = QuasiRefine(PRLog,N);

for i=1:nbParam
    % Get the value on log scale
    value = GetParam(PRLog,params{i});
    % Get the epsi on log scale
    position = FindParam(Sys,params(i)); %look for epsi index in ParamList
    position = PRLog.dim==position; %then for the line in the epsi array
    epsi = PRLog.epsi(position,:);
    
     % possibility 1 : define epsi, using the lower bound
%    inf = value - epsi; %lowest value of new intervalles on log scale
%    value = 10.^value; %value of new points on normal scale
%    inf = 10.^inf; %minimal value of the intervalles on normal scale
%    epsi = value - inf; %compute new epsi in normal scale
    
    %possiblity 2 : define epsi using the higher bound
%    sup = value + epsi; %highest value for new intervalles on log scale
%    value = 10.^value; %new points on normal scale
%    sup = 10.^sup; %highest value on normal scale
%    epsi = sup - value; %compute new epsi on normal scale
    
    %possibility 3 : define value using the lower and higher bounds
    %inf = value - epsi; %lowest new value on log scale
    inf = max( value-epsi , ranges(i,1) );
    %sup = value + epsi; %highest new value on log scale
    sup = min( value+epsi , ranges(i,2) );
    inf = 10.^inf; %lowest new value on normal scale
    sup = 10.^sup; %highest new value on normal scale
    epsi = (sup-inf)/2; %new epsi on normal scale
    value = inf + epsi; %compute new value on normal scale
    
    PRLog = SetParam(PRLog,params(i),value);
    PRLog.epsi(position,:) = epsi;
end

end
