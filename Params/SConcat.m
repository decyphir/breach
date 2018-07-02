function P = SConcat(P, P2)
%SCONCAT concatenates two parameter sets.
% 
% Synopsis: P = SConcat(P, P2)
% 
%  Tries to concat all compatible fields in P2 to those in P1. Basically
%  add points and trajectories of P2 in P1. If P2 contains property
%  parameters not in P1, there are added to P1. The value of these
%  parameters for the trajectories of P1 is set to 0. If a parameter is
%  uncertain for the parameter set Pi and not for Pj, they are uncertain in
%  the resulting parameter set. The epsi for these paramters for the
%  trajectories coming from Pj is 0.
% 
% Example (Lorentz84):
%   CreateSystem;
%   P1 = CreateParamSet(Sys,'a',[0.15,0.75],3);
%   P1 = SetParam(P1,'b',4);
%   P2 = CreateParamSet(Sys,'b',[3,9],3);
%   P2 = SetParam(P2,'a',0.25);
%   P = SConcat(P1,P2);  % contains 6 parameter vectors. Two are the same.
%   figure ; SplotBoxPts(P);
%   
%   P1 = ComputeTraj(Sys,P1,0:0.1:10);
%   P2 = ComputeTraj(Sys,P2,0:0.1:10);
%   P = SConcat(P1,P2);
%   numel(P.traj)  % It contains only 5 trajectories
% 
%See also Sselect
%

% Do nothing if one argument is empty
if isempty(P2)
    return;
end
if isempty(P)
    P = P2;
    return;
end

% check consistency between the to parameter sets
if(P.DimP~=P2.DimP)
    error('SConcat:DimP','The number of system parameters is not the same');
end
if(P.DimX~=P2.DimX)
    error('SConcat:DimX','The number of variables is not the same');
end
if ~isempty(setxor(P.ParamList(1:P.DimP),P2.ParamList(1:P2.DimP)))
    error('SConcat:ParamList','The system parameters are not the same');
end

% %%%%
% time_mult
% %%%%

if isfield(P2,'time_mult')
    if isfield(P,'time_mult')
        if(P.time_mult~=P2.time_mult)
            warning('SConcat:time_multField','time_mult field differs, keeping first one');
        end
    else
        P.time_mult = P2.time_mult;
    end
end


% %%%%
% props, props_names, props_values
% %%%%

if(~isfield(P,'props_names') || ~isfield(P2,'props_names'))
    try %#ok<TRYNC>
        P = rmfield(P,'props_names');
        P = rmfield(P,'props');
        P = rmfield(P,'props_values');
    end
else
    [P.props_names,iP,iP2] = intersect(P.props_names,P2.props_names,'stable');
    P.props = P.props(iP);
    P.props_values = [P.props_values(iP,:),P2.props_values(iP2,:)];
end

field_list = {'XS0', 'Xf', 'ExpaMax', 'XSf'};

for ii = 1:numel(field_list)
    if(isfield(P,field_list{ii}) && isfield(P2,field_list{ii}))
        if(numel(P.(field_list{ii}))==0)
            P.(field_list{ii}) = P2.(field_list{ii});
        else
            P.(field_list{ii}) = [P.(field_list{ii}) P2.(field_list{ii})];
        end
    end
end


% %%%%
% selected
% %%%%

if ~isfield(P,'selected')
    nParam = 0;
    if isfield(P,'pts')
        nParam = size(P.pts,2);
    end
    P.selected = false(1,nParam);
end
if ~isfield(P2,'selected')
    nParam = 0;
    if isfield(P2,'pts')
        nParam = size(P2.pts,2);
    end
    P2.selected = false(1,nParam);
end
P.selected = [P.selected,P2.selected];

if isempty(P2.pts)
    return;
end

if isempty(P.pts)
    P = P2;
    return;
end


% %%%%
% traj, traj_ref
% %%%%

if isfield(P,'traj') % Some check, in case of ...
    if(size(P.traj,1) > 1)
        P.traj = P.traj';
    end
    if ~isfield(P, 'traj_ref')
        P.traj_ref = 1:numel(P.traj);
    end
end


if isfield(P2,'traj')
    if(size(P2.traj,1) > 1)
        P2.traj = P2.traj';
    end
    if ~isfield(P2, 'traj_ref')
        P2.traj_ref = 1:numel(P2.traj);
    end
end

if(isfield(P,'traj') && isfield(P2,'traj'))
    num_traj_P = numel(P.traj);
    P2.traj_ref(P2.traj_ref~=0) = P2.traj_ref(P2.traj_ref~=0) + num_traj_P;
    
    % link param vector of P2 to traj of P
    for ii = 1:num_traj_P
        P2.traj_ref(ismember(P2.pts(1:P2.DimP,:)',P.traj{ii}.param,'rows')) = ii; % P.traj{ii}.param is a row vector
    end
    
    % copy P2.traj not in P.traj
    [traj_valid,~,i_unique] = unique(P2.traj_ref(P2.traj_ref>num_traj_P),'stable');
    i_unique = reshape(i_unique,1,[]);
    P.traj = [ P.traj , P2.traj(traj_valid-num_traj_P) ];
    
    % update P2.traj_ref (for traj not in P)
    P2.traj_ref(P2.traj_ref>num_traj_P) = i_unique+num_traj_P;
    
    % link param vector of P to P2 trajectories
    for ii=num_traj_P+1:numel(P.traj)
        P.traj_ref(ismember(P.pts(1:P.DimP,:)',P.traj{ii}.param,'rows')) = ii;
    end
    
    % copy P2.traj_ref in P.traj_ref
    P.traj_ref = [P.traj_ref, P2.traj_ref];
elseif isfield(P,'traj')
    % link param vector of P2 to traj of  P
    P2.traj_ref = zeros(1,size(P2.pts,2));
    for ii = 1:numel(P.traj)
        P2.traj_ref(ismember(P2.pts(1:P2.DimP,:)',P.traj{ii}.param,'rows')) = ii; % P.traj{ii}.param is a row vector
    end
    % copy P2.traj_ref in P.traj_ref
    P.traj_ref = [P.traj_ref, P2.traj_ref];
elseif isfield(P2,'traj')
    % link param vector of P to traj of  P2
    P.traj_ref = zeros(1,size(P.pts,2));
    for ii = 1:numel(P2.traj)
        P.traj_ref(ismember(P.pts(1:P.DimP,:)',P2.traj{ii}.param,'rows')) = ii; % P2.traj{ii}.param is a row vector
    end
    % copy traj
    P.traj = P2.traj;
    % copy P2.traj_ref
    P.traj_ref = [P.traj_ref, P2.traj_ref];
else
    P.traj_ref = zeros(1,size(P.pts,2)+size(P2.pts,2));
end


% %%%%
% pts
% %%%%


if(isfield(P,'pts') && isfield(P2,'pts'))  % case where P.pts or P2.pts are empty already considered
    % for pts in P, we set params of P2\P to last values in P2
    newParams=setdiff(P2.ParamList,P.ParamList,'stable');
    if ~isempty(newParams)
        val = GetParam(P2, newParams);
        P=SetParam(P,newParams,val(end,:));
    end
    
    % for pts in P2, we set the params of P\P2 to last value in P
    newParams=setdiff(P.ParamList,P2.ParamList,'stable');
    if ~isempty(newParams)
        val = GetParam(P, newParams);
        P2=SetParam(P2,newParams,val(:,end));
    end
    
    % we add P2.pts to P.pts
    P.pts = [ P.pts , P2.pts(FindParam(P2,P.ParamList),:) ];
end


% %%%%
% epsi, dim
% %%%%

try % keeps backward compat. to some extent
if(isfield(P,'epsi') && isfield(P2,'epsi'))
    
    % Add new uncertain parameter in P
    name_new = setdiff(P2.ParamList(P2.dim),P.ParamList(P.dim));
    P.dim = [P.dim, FindParam(P,name_new)];
    P.epsi = [P.epsi;zeros(numel(name_new),size(P.epsi,2))];
    
    % Add new uncertain parameters in P2
    name_new = setdiff(P.ParamList(P.dim),P2.ParamList(P2.dim));
    P2.dim = [P2.dim, FindParam(P2,name_new)];
    P2.epsi = [P2.epsi;zeros(numel(name_new),size(P2.epsi,2))];
    
    % get index of P2.dim in P.dim
    [~,~,iP2] = intersect(P.ParamList(P.dim), P2.ParamList(P2.dim),'stable'); % P.ParamList(P.dim) and P2.ParamList(P2.dim) are the same
    
    % Copy P2.epsi in P.epsi
    P.epsi = [P.epsi, P2.epsi(iP2,:)];
end

catch % don't sweat it
    P.epsi = zeros(numel(P.dim), size(P.pts, 2)); 
end

% %%%%
% traj_to_compute
% %%%%

[~,P.traj_to_compute] = unique(P.pts(1:P.DimP,:)','rows','first');
if isfield(P,'traj') % optimizing test
    P.traj_to_compute = setdiff(P.traj_to_compute,find(P.traj_ref~=0)); % don't keep those already computed
end
P.traj_to_compute = sort(reshape(P.traj_to_compute,1,[]));


% %%%%
% init_fun
% %%%%

if(isfield(P2,'init_fun') && ~isfield(P,'init_fun'))
    P.init_fun = P2.init_fun;
end

end
