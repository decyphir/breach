function Pf = ComputeTrajSensi(Sys, P0, tspan, i_params)
%COMPUTETRAJSENSI computes trajectories with corresponding sensitivities
% issued from points in P0 on the time interval tspan. The function
% recomputes all trajectories, even if they were already computed.
% 
% Synopsis:  Pf = ComputeTrajSensi(Sys, P0, tspan[, i_params])
% 
% Inputs:
%   -  Sys      : System (needs to be compiled)
%   -  P0       : Initial parameter set
%   -  tspan    : interval of the form [t0, tf], t0:dt:tf, etc
%   -  i_params : Names or indexes of parameter sensitivities to compute,
%                 if absent uses uncertain parameters in P (aka P.dim). Not
%                 valid indexes or names will not be considered.
% 
% Output:
%   -  Pf : Parameter set augmented with the field traj containing
%           computed trajectories with sensitivities
% 
% Example (Lorentz84):
%  CreateSystem;
%  P = CreateParamSet(Sys,'x0',[-1,1],2);
%  P = ComputeTrajSensi(Sys,P,0:0.1:10)
%  %TO FINISH
% 
%See also ComputeTraj SpropSensi PPhiSensiLocal SplotSensiBar
%

%%%%%%%%%%%%%%%%%%
% Manage inputs

if isfield(Sys,'type')
    if strcmp(Sys.type,'traces')
        Pf = P0;
        return ;
    end
else
    Sys.type = 'Breach';
end

% manage i_params
if ~exist('i_params','var')
    i_params = [];
elseif(iscell(i_params) || ischar(i_params))
    i_params = FindParam(P0,i_params);
end
i_params = i_params(i_params<=size(P0.pts,1));
i_params = i_params(i_params>0);

% do initialization if such function exists
if ~isempty(i_params)
    org_dims = P0.dim; % we save original dim and epsi
    org_epsi = P0.epsi;
    P0 = SAddUncertainParam(P0,i_params);
    P0 = SDelUncertainParam(P0,org_dims,i_params);
end
if isfield(Sys, 'init_fun')
    P0 = Sys.init_fun(P0);
end
if isfield(P0, 'init_fun')
    P0 = P0.init_fun(P0);
end

[~,traj_to_compute] = unique(P0.pts(1:P0.DimP,:)','rows','first');
if(numel(traj_to_compute)~=size(P0.pts,2)) % we got duplicate
    traj_to_compute = sort(traj_to_compute);
    Ptmp = Sselect(P0,traj_to_compute);
    Ptmp = ComputeTrajSensi(Sys,Ptmp,tspan,i_params);
    
    Pf = P0;
    Pf.traj = Ptmp.traj; % copy the computed fields
    Pf.Xf = Ptmp.Xf;
    Pf.XSf = Ptmp.XSf;
    Pf.ExpaMax = Ptmp.ExpaMax;
    [~,Pf.traj_ref] = ismember(Pf.pts(1:Pf.DimP,:)',Ptmp.pts(1:Ptmp.DimP,:)','rows'); % link all param vect to computed traj and sensi
    Pf.traj_ref = reshape(Pf.traj_ref,1,[]); % set traj_ref in a line shape
    Pf.traj_to_compute = [];
    return;
end

% manage tspan
if iscell(tspan)
    if(numel(tspan)==2)
        T = [tspan{1} tspan{2} tspan{2}];
    else
        T = cell2mat(tspan);
    end
else
    T = tspan;
end

%%%%%%%%%%%%%%%%%%
% Do the computation

InitSensi(Sys,P0);

if(isfield(P0,'XS0') && isempty(P0.XS0))
    P0 = rmfield(P0,'XS0');
end

if ~isfield(P0,'XS0') % XS0 describes the sensitivity of each variable to each uncertain variable and param
    dims = sort(P0.dim);
    Ns = numel(dims);
    N = P0.DimX;
    ix0 = dims(dims<=N); % indexes in pts of uncertain initial conditions
    
    xS0 = zeros(N,Ns);
    
    for ii=1:numel(ix0); % for each uncertain initial condition
        xS0(dims(ii),ii) = 1; % d x_ii(0)/d x_ii = 1
    end
    
    xS0 = reshape(xS0,N*Ns,1); % we stack columns of yS0
    
    P0.XS0 = repmat(xS0,[1 size(P0.pts,2)]); % duplicate initial sensitivity for each parameter vector
    
end

Pf = cvm(93,P0,T); % fill the fields traj, XS (?) and Xf (?)

CVodeFree();

if ~isempty(i_params)
    Pf = SAddUncertainParam(Pf,org_dims);
    Pf = SDelUncertainParam(Pf,i_params,org_dims);
    Pf.epsi = org_epsi;
end

Pf.traj_to_compute = [];
Pf.traj_ref = 1:size(P0.pts,2);

end
