function ParamValues = GetParam(S,ParamList)
% GETPARAM answers a matrix containing the values of parameters in the
% parameter set S. Each line contains the values of one parameter. The
% function returns an empty matrix if a parameter name or a parameter
% indice is not valid. S can be a system or a parameter set. If S is a
% system, GetParam answers the value provided when the system had been
% created.
%
% Synopsis: ParamValues = GetParam(S, ParamList)
%
% Input :
%   - S         : a system or a parameter set
%   - ParamList : either a cell matrix containing the name of the
%                 parameters to get, or an array containing the indexes of
%                 the parameters, or a string containing the name of the
%                 parameter.
%
% Example (for Lorenz84 system):
%
%    CreateSystem;
%    P = CreateParamSet(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    val = GetParam(Pr, 'a');
%    val = GetParam(Pr, 'F'); % get an array 1x9 filled with the value of F
%    val = GetParam(Pr, 'blah'); % return an empty array
%    val = GetParam(Pr, {'a', 'b'});
%    val = GetParam(Pr, {'b', 'a'});
%
%See also SetParam, GetEpsi

    if (isfield(S,'pts'))
        p = S.pts;
    else
        p = S.p;
    end

    % case 1 : the ParamList argument is an array char of one parameter
    % name
    if ischar(ParamList)
        ind = FindParam(S,ParamList);
        if(ind<=size(p,1)) %check parameter existence
            ParamValues = p(ind,:);
        else
            ParamValues = [];
        end
        return
    end

    % case 2 : the ParamList parameter is a list of integer
    if isnumeric(ParamList)
        %check for an error in the list
        if ~isempty(find(ParamList>size(p,1),1)) || ...
                                    ~isempty(find(ParamList<=0,1))
            ParamValues = [];
            return
        end
        %initialize ParamValues
        ParamValues = zeros(numel(ParamList),size(p,2));
        %copy the values
        for i = 1:numel(ParamList)
            ParamValues(i,:) = p(ParamList(i),:);
        end
        return
    end

    % case 3 : the ParamList parameter is a list of parameter name
    inds = FindParam(S,ParamList);
    %initialization
    if ~isempty(find(inds>size(p,1),1))
        ParamValues = [];
        return;
    else
        ParamValues = zeros(numel(ParamList),size(p,2));
    end
    for i = 1:numel(ParamList)
        ParamValues(i,:) = p(inds(i),:);
    end

end
