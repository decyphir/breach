function ParamValues = GetParam(S,ParamList)
%GETPARAM returns a matrix containing the values of parameters in S.
% Each line contains the values of one parameter. The function returns an
% empty matrix if a parameter name or a parameter indice is not valid. S
% can be a system or a parameter set. If S is a system, GetParam returns
% the value provided when the system has been created.
%
% Synopsis: ParamValues = GetParam(S, ParamList)
%
% Inputs:
%   - S         : a system or a parameter set
%   - ParamList : either a cell array containing the names of the
%                 parameters to get, or an array containing the indexes of
%                 the parameters, or a string containing the name of the
%                 parameter.
%
% Output:
%  - ParamValues : an array of size numel(ParamList) x size(S.p,2), except
%                  if ParamList is a char array in which case, its size is
%                  1 x size(S.p,2). Consider S.pts if S is a parameter set
%                  instead of a system.
%
% Example (for Lorenz84 system):
%
%    CreateSystem;
%    P = CreateParamSet(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    val = GetParam(Pr, 'a')
%    val = GetParam(Pr, 'F') % get an array 1x9 filled with the value of F
%    val = GetParam(Pr, 'blah') % return an empty array
%    val = GetParam(Pr, {'a', 'b'})
%    val = GetParam(Pr, {'b', 'a'})
%    val = GetParam(Pr, {'a', 'blah'}) % return an empty array
%
%See also SetParam, GetEpsi

if isfield(S,'pts')
    p = S.pts;
else
    p = S.p;
end

ParamValues = []; % default value

% ParamList is a list of integer
if isnumeric(ParamList)
    %check for an error in the list
    if all(ParamList<=size(p,1)) && all(ParamList>0)
        ParamValues = p(ParamList,:);
    end
    
% ParamList is a parameter name or a cell array of parameter names
elseif(ischar(ParamList) || iscell(ParamList))    
    inds = FindParam(S,ParamList);
    %check for an error in the list
    if all(inds<=size(p,1))
        ParamValues = p(inds,:);
    end
end

end
