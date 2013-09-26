function Pf = ComputeTrajSensi(Sys, P, tspan, is)
%COMPUTETRAJSENSI computes trajectories with corresponding sensitivities
% issued from points in P0 on the time interval tspan.
%
% Synopsis:  Pf = ComputeTrajSensi(Sys, P0, tspan[, is])
%
% Inputs:
%
%    -  Sys   : System (needs to be compiled)
%    -  P0    : Initial parameter set
%    -  tspan : interval of the form [t0, tf], t0:dt:tf, etc
%    -  is    : Names or indexes of parameter sensitivities to compute, if
%               absent uses uncertain parameters in P (aka P.dim). Not
%               valid indexes or names will not be considered.
%
% Output:
%
%    -  Pf : Parameter set augmented with the field traj containing
%            computed trajectories with sensitivities
%
%  NEEDS UPDATE TO HANDLE TRAJ_REF AND PROPERTIES
%


if isfield(Sys, 'type')
    if strcmp(Sys.type,'traces')
        Pf = P;
        return ;
    end
end

% manage inputs
if iscell(tspan)
    if(numel(tspan)==2)
        T = [tspan{1} tspan{2} tspan{2}];
    else
        T = cell2mat(tspan);
    end
else
    T = tspan;
end

if ~exist('is','var')
    is = [];
elseif(iscell(is) || ischar(is))
    is = FindParam(is);
end
is = is(is<size(P.pts,2));
    

% checks for an initialization function
if isfield(Sys,'init_fun')
    P = Sys.init_fun(P);
end
if isfield(P,'init_fun')
    P = P.init_fun(P);
end

if ~isempty(is)
    org_dims = P.dim; % we save original dim and epsi
    org_epsi = P.epsi;
    P = SAddUncertainParam(P,is);
    P = SDelUncertainParam(P,org_dims,is);
end

InitSensi(Sys,P);

if(isfield(P,'XS0') && isempty(P.XS0))
    P = rmfield(P,'XS0');
end

if ~isfield(P,'XS0') % XS0 describes the sensitivity of each variable to each uncertain variable and param
    dims = P.dim;
    Ns = numel(dims);
    N = P.DimX;
    ix0 = dims(dims<=N); % indexes in pts of uncertain initial conditions
    
    yS0 = zeros(N,Ns);
    
    for ii=1:numel(ix0); % for each uncertain initial condition
        yS0(dims(ii),ii) = 1;
    end
    
    xS0 = reshape(yS0,N*Ns,1); % we stack columns of yS0
    
    P.XS0 = repmat(xS0,[1 size(P.pts,2)]); % duplicate initial sensitivity for each parameter vector
    
end

Pf = cvm(93,P,T); % fill the fields traj, XS (?) and Xf (?)
CVodeFree();

if ~isempty(is)
    Pf = SAddUncertainParam(Pf,org_dims);
    Pf = SDelUncertainParam(Pf,is,org_dims);
    Pf.epsi = org_epsi; %NM : not sure of that (wonder if cvm makes changes on epsi - I guess no)
end

end
