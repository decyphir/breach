function S = SConcat(S,S2)
% SCONCAT Concatenates two parameter sets.
%
% Synopsis: P = SConcat(P1, P2)
%
%  Tries to concat all compatible fields in P2 to those in P1. Basically
%  add points and trajectories of P2 in P1. If P2 contains property
%  parameters not in P1, there are added to P1. The value of these
%  parameters for the trajectories of P1 is set to 0. If a parameter is
%  uncertain for the parameter set Pi and not for Pj, they are uncertain in
%  the resulting parameter set. The epsi for these paramters for the
%  trajectories coming from Pj is 0.
%
% Example:
%   CreateSystem;
%   P1 = CreateParamSet(Sys,'a',[0.2,0.6],2);
%   P2 = CreateParamSet(Sys,'b',[3.5,7.5],2);
%   P = SConcat(P1,P2);
%   SplotBoxPts(P);
%


% check consistency between the to parameter sets
if(S.DimP~=S2.DimP)
    error('SConcat:DimP','The number of system parameters is not the same');
end

if ~isempty(setdiff(S.ParamList(1:S.DimP),S2.ParamList(1:S2.DimP)))
    error('SConcat:ParamList','The system parameters are not the same');
end


field_list_copy = {'props_names', 'props',  'time_mult'};

for i = 1:numel(field_list_copy)
    if isfield(S2, field_list_copy{i})
        S.(field_list_copy{i}) = S2.(field_list_copy{i});
    end
end

field_list = {'XS0', 'Xf', 'ExpaMax', 'XSf', 'props_values'};

for i = 1:numel(field_list)
    
    if ((isfield(S,field_list{i}))&&(isfield(S2,field_list{i})))
        if (numel(S.(field_list{i}))==0)
            S.(field_list{i}) = S2.(field_list{i});
        else
            S.(field_list{i}) = [S.(field_list{i}) S2.(field_list{i})];
        end
    end
end

% %%%%
% selected
% %%%%

if isfield(S2,'selected') && isfield(S,'pts')
    S2.selected = S2.selected+size(S.pts,2);
    if isfield(S,'selected')
        S.selected = [S.selected,S2.selected];
    else
        S.selected = S2.selected;
    end
end

% %%%%
% pts
% %%%%

if isfield(S,'pts') && isfield(S2,'pts')
    % for pts in S, we set to 0 the params of S2\S
    newParams=setdiff(S2.ParamList,S.ParamList);
    S=SetParam(S,newParams,zeros(1,numel(newParams)));
    
    % number of points in S
    nb_pts_S = size(S.pts,2);
    
    % for pts in S2, we set to 0 the params of S\S2
    S.pts = [S.pts,zeros(size(S.pts,1),size(S2.pts,2))];
    
    % we copy S2 into S
    for i=1:numel(S2.ParamList)
        idx = FindParam(S,S2.ParamList(i));
        S.pts(idx,nb_pts_S+1:end) = S2.pts(i,:);
    end
end

% %%%%
% epsi and dim
% %%%%

if isfield(S,'epsi') && isfield(S2,'epsi')
    % there are size(S2.pts,2) more trajectories
    S.epsi = [S.epsi,zeros(size(S.epsi,1),size(S2.pts,2))];
    
    for i=1:numel(S2.dim) % for each uncertain parameter in S2
        p = S2.ParamList(S2.dim(i));
        idx_p = find(strcmp(S.ParamList(S.dim),p),1);
        if isempty(idx_p) % if it is not an uncertain parameter in S
            S.dim = [S.dim,find(strcmp(S.ParamList,p),1)]; % it becomes one
            S.epsi = [S.epsi;zeros(1,size(S.epsi,2))]; % the epsi are null
            S.epsi(size(S.epsi,1),nb_pts_S+1:end) = S2.epsi(i,:); % exept for those in S2
        else
            S.epsi(idx_p,nb_pts_S+1:end) = S2.epsi(i,:);
        end
    end
end



if (isfield(S, 'traj')&&(isfield(S2,'traj')))
    
    %%% TEMPORARY FIX
    if size(S.traj,1) >1
        S.traj = S.traj';
    end
    
    if size(S2.traj,1) >1
        S2.traj = S2.traj';
    end
    
    if (~isfield(S, 'traj_ref'))
        S.traj_ref= 1:numel(S.traj);
    end
    
    if (~isfield(S2, 'traj_ref'))
        S2.traj_ref= 1:numel(S2.traj);
    end
    
    S.traj_ref = [S.traj_ref S2.traj_ref+numel(S.traj)];
    S.traj = [S.traj S2.traj];
    
end

