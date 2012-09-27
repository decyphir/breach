function ParamValues = GetParam(S,ParamList)
% GETPARAM get the values of parameters in a parameter set. The function
% return an empty set if a parameter name or a parameter indice is not
% valid.
%
% Synopsis: ParamValues = GetParam(P, ParamList)
%
% Example (for Lorenz84 system):
%
%    CreateSystem;
%    P = CreateSampling(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    val = GetParam(Pr, 'a' );
%    val = GetParam(Pr, 'b' );
%    val = GetParam(Pr, {'a', 'b'});
%    val = GetParam(Pr, {'b','a'});
%

    if (isfield(S,'pts'))
        p = S.pts;
    else
        p = S.p;
    end

    % case 1 : the ParamList parameter is an array char of one parameter
    % name
    if ischar(ParamList)
        ind = FindParam(S,ParamList);
        if(ind<=S.DimP) %check parameter existence
            ParamValues = p(ind,:);
        else
            ParamValues = [];
        end
        return
    end

    % case 2 : the ParamList parameter is a list of integer
    if isnumeric(ParamList)
        %check for an error in the list
        if ~isempty(find(ParamList>S.DimP,1)) || ...
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
    if ~isempty(find(inds>S.DimP,1))
        ParamValues = [];
        return;
    else
        ParamValues = zeros(numel(ParamList),size(p,2));
    end
    for i = 1:numel(ParamList)
        ParamValues(i,:) = p(inds(i),:);
    end

end
