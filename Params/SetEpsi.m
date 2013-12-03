function P = SetEpsi(P, ParamList, EpsiValues)
%SETEPSI sets the values of epsi in a parameter set.
% 
% Synopsis: P = SetEpsi(P, ParamList, EpsiValues)
% 
% Inputs:
%  - P         : the parameter set to modify ;
%  - ParamList : the list of parameter names or indexes for which the epsi
%                value is modified. If empty, nothing is done. If a
%                parameter of ParamList is not uncertain or doesn't exist,
%                then it is skipped ;
%  - EpsiValue : the values of the epsi. Its size is either
%                numel(ParamList) x size(P.pts,2) or numel(ParamList) x 1
% 
% Output:
%  - P : the new parameter set
% 
% Example (Lorenz84):
%    CreateSystem;
%    P = CreateParamSet(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    val = GetEpsi(Pr, 'a');
%    val10 = 10*val; 
%    Pr10 = SetEpsi(Pr, 'a', val10); % epsi for 'a' in Pr10 are ten
%                                     % times those in Pr
% Other example :
% 
%    CreateSystem;
%    P = CreateParamSet(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    Pr_2 = SetEpsi(Pr, {'a','G','blah'}, [0.4; 0.6 ; 2]);
%                     % epsilon for 'a'is set to 0.4 for the nine parameter
%                     % sets. Nothing has changed for 'G'.
%    
%See also GetEpsi CreateParamSet SetParam SAddUncertainParam
%SDelUncertainParam
%

if isempty(ParamList)
    return ;
end

if(size(EpsiValues,1)==1 && ~ischar(ParamList) && numel(ParamList)~=1)
    % in case when we have one value for many parameters
    EpsiValues = EpsiValues';
end

if(ischar(ParamList) || iscell(ParamList))
    ParamList = FindParam(P,ParamList);
end

[~,ind,ind_epsi] = intersect(P.dim,ParamList);
P.epsi(ind,:) = EpsiValues(ind_epsi,:);

end
