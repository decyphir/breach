function P = Preset_traj_ref(P)
% Preset_traj_ref reset traj_ref and traj_to_compute fields 

P.traj_ref = zeros(1,size(P.pts,2)); % initialise traj_ref

%% If traj exists, 
if isfield(P,'traj')&&~isempty(P.traj)
    param_traj = zeros(numel(P.traj{1}.param),numel(P.traj));
    for itraj = 1:numel(P.traj)
        param_traj(:,itraj)=P.traj{itraj}.param'; % param value for computed traj
    end
    [~, iu] = unique(param_traj', 'rows', 'stable');
    idx_new_traj = 1;
    for ii = 1:numel(iu) % for each existing traj
        same = all(bsxfun(@eq,P.pts(1:P.DimP,:),param_traj(:,iu(ii))),1);
        if any(same) % if it is equal to a parameter vector
            if idx_new_traj~=ii 
                P.traj{idx_new_traj} = P.traj{iu(ii)}; % then keep the traj (ok because idx_new_traj<=ii)
            end
            P.traj_ref(same) = idx_new_traj;
            idx_new_traj = idx_new_traj + 1;
        end
    end
    P.traj = P.traj(1:idx_new_traj-1); % remove all other traj
    num_traj = numel(P.traj);

    if isempty(P.traj)
        P = rmfield(P,'traj');
    end
else
    num_traj=0;
end

pts_indices_to_compute = find(~logical(P.traj_ref));
if isempty(pts_indices_to_compute)
    P.traj_to_compute = [];
else
    
    [~, ia, ib] = unique(P.pts(1:P.DimP,pts_indices_to_compute)','rows','stable');
    P.traj_to_compute = pts_indices_to_compute(ia);
    traj_ref_to_compute = num_traj+(1:numel(P.traj_to_compute));
    
    P.traj_ref(pts_indices_to_compute) = ...   % this is where traj_ref is 0
            traj_ref_to_compute(ib);    

end