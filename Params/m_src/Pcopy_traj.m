function P = Pcopy_traj(P, P0)
%Pcopy_traj get trajectories from P0 in P, update traj_ref and traj_to_compute 
%

n = numel(P0.dim);
P.dim = P0.dim;

if isfield(P0,'traj_plot_opt')
    P.traj_plot_opt = P0.traj_plot_opt;
end

if isfield(P0,'X0plot_opt')
    P.X0plot_opt = P0.X0plot_opt;
end

if isfield(P0,'plot_proj')
    P.plot_proj = P0.plot_proj;
end

if isfield(P0,'ParamList')
    P.ParamList = P0.ParamList;
end

if isfield(P0,'selected')
    P.selected = zeros(1, size(P.pts,2));
end

if isfield(P0,'time_mult')
    P.time_mult = P0.time_mult;
end

P.DimX = P0.DimX;
P.DimP = P0.DimP;

%  Checks for pts sharing the sames systems parameters
P.traj_ref = zeros(1,size(P.pts,2));

% check if some already computed traj are valid for generated param set
[Ppts,~,i_pts] = unique(P.pts(1:P.DimP,:)','rows','stable'); % got unique param vector for P
[~,iP,iP0] = intersect(Ppts,P0.pts(1:P.DimP,logical(P0.traj_ref))','rows','stable'); % get param vector common between computed P0 and (unique) P

for ii = 1:numel(iP) % for all param vector common between P and those computed in P0
    P.traj_ref(i_pts==iP(ii)) = ii; % affect the ith traj to all param vector equals to the one in P0
    P.traj{ii} = P0.traj{P0.traj_ref(iP0(ii))}; % and fill the iith traj in P
    P.Xf(1:P.DimX,ii) = P0.Xf(1:P0.DimX,P0.traj_ref(iP0(ii)));
end

% Then define which are remaining trajectories to compute
[~,P.traj_to_compute] = unique(P.pts(1:P.DimP,:)','rows','first');
if isfield(P,'traj')
    P.traj_to_compute = setdiff(P.traj_to_compute,find(P.traj_ref~=0)); % don't keep those already computed
end
P.traj_to_compute = sort(reshape(P.traj_to_compute,1,[])); % set it in a line shape

end
