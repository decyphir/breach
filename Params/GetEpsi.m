function [EpsiValues, EpsiNames] = GetEpsi(P,ParamList)
%GETEPSI gets the epsi of parameters in a parameter set. The function does
% not provide a value for parameter name or indice not valid or not
% uncertain.
% 
% Synopsis: [EpsiValues, EpsiNames] = GetEpsi(P, ParamList)
% 
% Inputs:
%  - P         : the parameter set. It may contain many parameter vectors.
%  - ParamList : the list of parameter for which we get the epsilon. If it
%                is empty, nothing is done.
%
% Outputs:
%  - EpsiValues : an array of size numel(EpsiNames) x size(P.pts,2)
%                 containing the values of the epsi for the uncertain
%                 parameters. Each line contains the epsi for one
%                 parameter.
%  - EpsiNames  : cell array containing the name of uncertain parameters
%                 for which the epsilon has been provided in EpsiValue. It
%                 is empty if ParamList is empty or contains only not valid
%                 or fixed parameters.
% 
% Example (Lorenz84):
%   CreateSystem;
%   P = CreateParamSet(Sys, {'a','b'}, [0,9;0,5]);
%   P = Refine(P, 3);
%   val = GetEpsi(P, 'a') % epsi value of a is always the same (ie: 1.5)
%   [val,names] = GetEpsi(P, {'F','b','blah'}) % names is 'b', val contains
%                                             % its epsi value
% 
%See also GetParam SetEpsi SAddUncertainParam SDelUncertainParam
%

EpsiValues = [];
EpsiNames = {};

if(ischar(ParamList) || iscell(ParamList))
    ParamList = FindParam(P,ParamList);
end

if isempty(ParamList)
    return; % here, we manage empty array as well as empty cell array
end

%keep only valid parameters
[valid,idx_epsi] = ismember(ParamList,P.dim);
ParamList = ParamList(valid);
idx_epsi = idx_epsi(valid); % indexes in P.epsi
%set EpsiNames and EpsiValues
EpsiNames = P.ParamList(ParamList);
EpsiValues = P.epsi(idx_epsi,:);

end
