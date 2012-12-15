function EpsiValues = GetEpsi(P,ParamList)
% GETEPSI get the epsi of parameters in a parameter set. The function
% return an empty set if a parameter name or a parameter indice is not
% valid.
%
% Synopsis: EpsiValues = GetEpsi(P, ParamList)
%
%See also GetParam

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


    % case 1 : the ParamList argument is an array char of one parameter
    % name
    if ischar(ParamList)
        ind = FindParam(P,ParamList);
        ind = P.dim(P.dim==ind);
        if ~isempty(ind) %check epsi existence
            EpsiValues = P.epsi(ind,:);
        else
            EpsiValues = [];
        end
        return
    end

    % case 2 : the ParamList parameter is a list of integer
    if isnumeric(ParamList)
        %check if we ask for unexisting or not unknown parameters
        if ~isempty(setdiff(ParamList,P.dim))
            EpsiValues = [];
            return
        end
        %initialize ParamValues
        EpsiValues = zeros(numel(ParamList),size(P.pts,2));
        %copy the values
        for i = 1:numel(ParamList)
            ind = P.dim(P.dim==ParamList(i));
            EpsiValues(i,:) = P.epsi(ind,:);
        end
        return
    end

    % case 3 : the ParamList parameter is a list of parameter name
    inds = FindParam(P,ParamList);
    %initialization
    if ~isempty(setdiff(inds,P.dim))
        EpsiValues = [];
        return;
    else
        EpsiValues = zeros(numel(ParamList),size(P.pts,2));
    end
    for i = 1:numel(ParamList)
        ind = P.dim(P.dim==inds(i));
        EpsiValues(i,:) = P.epsi(ind,:);
    end

end
