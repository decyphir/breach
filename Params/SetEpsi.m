function P = SetEpsi(P, ParamList, EpsiValues)
% SETPARAM Sets the values of epsi in a parameter set.
% 
% Synopsis: P = SetParam(P, ParamList, EpsiValues)
% 
% Input:
%  - P          : the parameter set to modify ;
%  - ParamList  : the list of parameters for which the epsi value is
%                 modified. If empty, nothing is done ;
%  - EpsiValue  : the values of the epsi. Its size is either
%                 ( numel(ParamList) , size(P.pts,2) ) or
%                 ( numel(ParamList) , 1 )
% 
% Example (for Lorenz84 system):
% 
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
%See also GetEpsi CreateParamSet SetParam
%  

if isempty(ParamList)
    return ;
end

if(size(EpsiValues,1)==1 && ~ischar(ParamList) && numel(ParamList)~=1)
    % in case when we have one value for many parameters
    EpsiValues = EpsiValues';
end

if(ischar(ParamList) || iscell(ParamList))
    ind = FindParam(P,ParamList);
    [~,ind,ind_epsi] = intersect(P.dim,ind);
    P.epsi(ind,:) = EpsiValues(ind_epsi,:);
    return
% 
% elseif isnumeric(ParamList)
%     for i = 1:numel(ParamList)
%         if ParamList(i) > size(P.pts, 1)
%             error('SetParam:index','Index out of range')
%         end
%         P.pts(ParamList(i),:) = EpsiValues(i,:);
%     end
%     return
end

end
