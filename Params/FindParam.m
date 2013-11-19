function index =  FindParam(S, param)
% FINDPARAM Finds the indices of parameters given by their name for a given
% system or parameter set.
%
% Syntax: index = FindParam(Sys, param)
%     or  index = FindParam(P, param)
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
%
% Example (Lorentz84):
%
%  CreateSystem;
%  idx = FindParam(Sys, 'a')
%  numel(Sys.ParamList)  % 7 elements
%  idxs = FindParam(Sys, {'a','blah','F','bou'}) % 4  8  6  9 -> 'blah' and 'bou' are new params
%  idxs = FindParam(Sys, {'F','F'})
%  P = CreateParamSet(Sys,'a',[0,2]);
%  idxs = FindParam(P, {'bou','bou'})
%

% check inputs
if ~isfield(S,'ParamList')
    error('FindParam:noParamList','No parameter list ...');
end

if isempty(param)
    index = [];
    return;
elseif ~iscell(param)
    param = {param};
end

%NM: This implementation is very optimized: I have try to use ismember or
% intersect, it is 10 to 100 times slower.

index = zeros(1, numel(param));

newp = S.DimP;
if isfield(S, 'pts') % S.pts can be > S.DimP (due to property parameters)
    newp = size(S.pts,1);
end

for ii = 1:numel(param)
    same = strcmp(S.ParamList,param{ii});
    if any(same)
        index(ii) = find(same,1);
    else
        newp = newp+1;
        index(ii) = newp;
        S.ParamList = [S.ParamList, param(ii)];
    end
end

end
