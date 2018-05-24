function [index, type] =  FindParam(S, param)
% FINDPARAM (Legacy) Finds the indices of parameters given by their name for a given
% system or parameter set.
%
% Syntax: [index, status] = FindParam(Sys, param)
%     or  [index, status] = FindParam(P, param)
%
% Inputs:
%  - S     : is the system or a parameter set
%  - param : it can be a string or a cell of string when looking for more
%            than one parameter. If it contains many time the same
%            parameter name, then this index is answered as many
%            time than the parameter name is in param.
%
% Outputs:
%  - index  is an array of size 1 x numel(param) containing the index of
%           param in S.ParamList. It is greater than S.DimP for
%           parameter(s) not found. (for unfound parameters, the returned
%           indexes is compact, different for each parameter and the
%           smaller index is equal to the higher known parameter index + 1)
%  - status array of same size as index, contains 0 if not found, 1 if
%           signal, 2 if param, 3 if other
%

index = [];
type = [];    

% check inputs
if ~isfield(S,'ParamList')
    error('FindParam:noParamList','No parameter list ...');
end

sizep0 = S.DimP;
if isfield(S, 'pts') % S.pts can be > S.DimP (due to property parameters)
    sizep0 = size(S.pts,1);
end
sizep = sizep0;

if isempty(param)
    index = [];
    return;
elseif isnumeric(param)
    index = param;
    CheckIdx();
    return;
elseif ~iscell(param)
    param = {param};
end

index = zeros(1, numel(param));

for ii = 1:numel(param)
    same = strcmp(S.ParamList,param{ii});
    if any(same)
        index(ii) = find(same,1);
    else
        sizep = sizep+1;
        index(ii) = sizep;
    end
end

CheckIdx();

    function   CheckIdx()    
        type = 3*(index <= sizep0);
        type(index <= S.DimP) = 2;
        type(index<= S.DimX) = 1;
    end

end



