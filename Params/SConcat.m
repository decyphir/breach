function P = SConcat(P, P2, fast)
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
% init_fun
% %%%%

if(isfield(P2,'init_fun') && ~isfield(P,'init_fun'))
    P.init_fun = P2.init_fun;
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

field_list = {'Xf'};

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
% pts
% %%%%


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


% %%%%
% epsi, dim
% %%%%

if fast
     P.epsi = zeros(numel(P.dim), size(P.pts, 2)); 
else
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
end
% %%%%
% traj, traj_ref, traj_to_compute
% %%%%

if ~exist('fast', 'var')
    fast =false;
end


if(isfield(P,'traj') && isfield(P2,'traj'))

    if fast
        P.traj_ref = [P.traj_ref P2.traj_ref+numel(P.traj)]; % should be good enough, most common case is appending one trace
    else
       P = Preset_traj_ref(P);
    end
    P.traj = [P.traj P2.traj];
    
elseif isfield(P,'traj')
    P = Preset_traj_ref(P);    
elseif isfield(P2,'traj')
    % copy traj
    P.traj = P2.traj;  
    P = Preset_traj_ref(P);
else
    P = Preset_traj_ref(P);
end

end
