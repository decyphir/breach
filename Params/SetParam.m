function P = SetParam(P, ParamList, ParamValues)
%SETPARAM sets the values of parameters in a parameter set. Note that if
% the parameter is not present in P, it is created and appended as a fixed
% parameter. The epsi are neither modified nor created.
% 
% Synopsis: P = SetParam(P, ParamList, ParamValues)
% 
% Inputs:
%  - P          : the parameter set to modify;
%  - ParamList  : a cell array containting the list of parameter names for
%                 which the value is changed or created. If empty, nothing
%                 is done;
%  - ParamValue : matrix containing the values of the parameters. Its size
%                 is either ( numel(ParamList) , size(P.pts,2) ) or
%                 ( numel(ParamList) , 1 ), in which case, the same value
%                 is set to each parameter.
% 
% Output:
%  - P : the modified parameter set
% 
% Example (for Lorenz84 system):
%   CreateSystem
%   P = CreateParamSet(Sys, {'a', 'b'}, [0 9; 0 5]);
%   Pr = Refine(P, 3);
%   val = GetParam(Pr, 'a')
%   val10 = 10*val; 
%   Pr10 = SetParam(Pr, 'a', val10); % values for 'a' in Pr10 are ten
%                                    % times those in Pr
%   GetParam(Pr10,'a') % 10 times higher
%
%   Pr_2 = SetParam(Pr, {'F','G'}, [0.4; 0.6]); % values for 'F' (resp.
%                    % 'G') equals to 0.4 (resp 0.6) in the nine points.
%                    % Values of 'a' and 'b' are equals in Pr and Pr_2.
%   GetParam(Pr_2,{'a','F'})
%   
%   Pr_3 = SetParam(Pr, 'a', 5); % a is equal to 5 in the nine points
%   GetParam(Pr_3,'a')
%   
%   Pr = ComputeTraj(Sys, Pr, 0:0.1:10);
%   Pr.traj_ref              % should be [1:9]
%   Pr_4 = SetParam(Pr,'a',1.5);
%   Pr_4.traj_ref            % should be [1,1,1,2,2,2,3,3,3]
%   Pr_4.traj_to_compute     % should be empty
%   
%   Pr_5 = SetParam(Pr,'a',2);
%   Pr_5.traj_ref            % should be all zero
%   Pr_5.traj_to_compute     % should be [1,4,7]
% 
%See also GetParam CreateParamSet SetEpsi
%

%%%%%%
% manage input

if isempty(ParamList)
    return ;
elseif ischar(ParamList)
    ParamList = {ParamList};
end

if(size(ParamValues,1) == 1 && numel(ParamList)~=1)
    ParamValues = ParamValues'; 
end

if(numel(ParamList) ~= size(ParamValues,1) && size(ParamValues,1)~=1)
    error('SetParam:wrongElementNumber','The number of rows of ParamValues is not the same as the number of ParamList');
end

if isfield(P,'pts')        % if P is a parameter set
    pts = 'pts';
else
    pts = 'p';
end
if ~isfield(P,'traj_ref')
    P.traj_ref = zeros(1,size(P.(pts),2));
    if isfield(P,'traj') % there is a structural error in P: try to fix
        m = min(size(P.(pts),2),numel(P.traj));
        P.traj_ref(1:m) = 1:m;
    end
end

%%%%%%
% set params
if iscell(ParamList)
    inds = FindParam(P,ParamList);
    new_params = ParamList(inds>size(P.(pts),1));
    P.ParamList = [P.ParamList new_params];
    
    if (size(P.(pts), 2)==1 && size(ParamValues,2)>1)
        P.(pts) = repmat(P.(pts), [1 size(ParamValues,2)]);
        P.(pts)(inds,:) = ParamValues;               
        if isfield(P, 'epsi') % I guess we'll want to check other fields too,e.g., selected etc
            P.epsi = repmat(P.epsi, [1 size(ParamValues,2)]);
        end
    else
        for ii = 1:numel(inds) % doing a loop is faster than using repmat on ParamValues or using bsxfun
            P.(pts)(inds(ii),:) = ParamValues(ii,:);
        end
    end

elseif isnumeric(ParamList)
    if any(ParamList>size(P.(pts), 1))
        warning('SetParam:index','Indexes out of range have been skipped')
        ParamValues = ParamValues(ParamList<=size(P.(pts), 1),:);
        ParamList = ParamList(ParamList<=size(P.(pts), 1));
    end
    
    for ii = 1:numel(ParamList) % doing a loop is faster than using repmat on ParamValues
        P.(pts)(ParamList(ii),:) = ParamValues(ii,:);
    end
end


%%%%%%
% manage traj_ref and traj_to_compute
P.traj_ref = zeros(1,size(P.(pts),2)); % initialise traj_ref
if isfield(P,'traj')
    param_traj = cat(1,P.traj(:).param)'; % param value for computed traj
    idx_new_traj = 1;
    for ii = 1:size(param_traj,2) % for each existing traj
        same = all(bsxfun(@eq,P.(pts)(1:P.DimP,:),param_traj(:,ii)),1);
        if any(same) % if it is equal to a parameter vector
            P.traj(idx_new_traj) = P.traj(ii); % then keep the traj (ok because idx_new_traj<=ii)
            P.traj_ref(same) = idx_new_traj;
            idx_new_traj = idx_new_traj + 1;
        end
    end
    P.traj = P.traj(1:idx_new_traj-1); % remove all other traj
    if isempty(P.traj)
        P = rmfield(P,'traj');
    end
end

[~,P.traj_to_compute] = unique(P.(pts)(1:P.DimP,:)','rows','first'); % we define traj_to_compute
if isfield(P,'traj')
    P.traj_to_compute = setdiff(P.traj_to_compute,find(P.traj_ref~=0)); % don't keep those already computed
end
P.traj_to_compute = sort(reshape(P.traj_to_compute,1,[])); % we set P.traj_to_compute in a line shape

end
