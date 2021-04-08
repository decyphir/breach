function P = TestronRefine(P, nb, varargin)

dim_num = numel(P.dim);
old_epsi = P.epsi;
P.epsi = repmat(P.epsi,[1 nb]);

r = rand(dim_num, nb);

width = 2 * old_epsi;
mini = P.pts(P.dim,:) - old_epsi;

P.pts = repmat(P.pts,[1 nb]);
P.pts(P.dim,:) = repmat(width,[1 nb]).*r+repmat(mini,[1 nb]);

P = Preset_traj_ref(P);


end