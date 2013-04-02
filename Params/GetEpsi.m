function [EpsiValues, EpsiNames] = GetEpsi(P,ParamList)
% GETEPSI get the epsi of parameters in a parameter set. Each line contains
% the epsi of one parameter. The function return an empty set if a
% parameter name or a parameter indice is not valid.
%
% Synopsis: [EpsiValues, EpsiNames] = GetEpsi(P, ParamList)
%
% Input:
%  - P         : the parameter set
%  - ParamList : the list of parameter for which we get the epsilon. if it
%                is empty, nothing is done.
%
% Output:
%  - EpsiValues : the values of the epsilon for the uncertain parameters
%  - EpsiNames  : the name of uncertain parameters for which the epsilon
%                 has been provided in EpsiValue
%
% Example (for Lorenz84 system):
%
%    CreateSystem;
%    P = CreateParamSet(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    val = GetEpsi(Pr, 'a'); % here the epsi value of a is always the same
%    [val,epsi] = GetEpsi(Pr, {'F','b','blah'}); % epsi is 'b', and val
%                            % contains its epsilon value
%
%See also GetParam SetEpsi SAddUncertainParam SDelUncertainParam
%

if isempty(ParamList)
    EpsiValues = [];
    EpsiNames = {};
end

% case 1 : the ParamList argument is an array char of one parameter
% name
if ischar(ParamList)
    ind = FindParam(P,ParamList);
    ind = P.dim==ind;
    if ~isempty(ind) %check epsi existence
        EpsiValues = P.epsi(ind,:);
        EpsiNames = {ParamList};
    else
        EpsiValues = [];
        EpsiNames = {};
    end
    return

% case 2 : the ParamList parameter is a list of integer
elseif isnumeric(ParamList)
    %check if we ask for unexisting or not unknown parameters
    if ~isempty(setdiff(ParamList,P.dim))
        EpsiValues = [];
        return
    end
    %initialize EpsiValues
    EpsiValues = zeros(numel(ParamList),size(P.pts,2));
    %copy the values
    for i = 1:numel(ParamList)
        ind = P.dim==ParamList(i);
        EpsiValues(i,:) = P.epsi(ind,:);
    end
    return

% case 3 : the ParamList parameter is a list of parameter name
elseif iscell(ParamList)
    inds = FindParam(P,ParamList);
    [inds,~,i_dim] = intersect(inds,P.dim,'stable');
    %initialization
    if isempty(inds)
        EpsiValues = [];
        EpsiNames = {};
        return;
    end
    EpsiValues = P.epsi(i_dim,:);
    EpsiNames = P.ParamList(inds);
end

end
