function Pf = ComputeTrajSensi(Sys, P, tspan, is)
%COMPUTETRAJSENSI compute trajectories with sensitivities
%
% Synopsis:  Pf = ComputeTrajSensi(Sys, P0, tspan [ , is])
% 
% Compute trajectories with corresponding sensitivities issued from points
% in P0 on the time interval tspan.
%
% Inputs:
%
%    -  Sys      System (needs to be compiled)
%    -  P0       Initial parameter set
%    -  tspan    interval of the form [t0, tf], t0:dt:tf, etc
%    -  is       Parameter sensitivities to compute, if absent uses
%                uncertain parameters in P (aka P.dim)
%
% Outputs:
%
%    -  Pf       Parameter set augmented with the field traj
%                containing computed trajectories with sensitivities
%
%  NEEDS UPDATE TO HANDLE TRAJ_REF AND PROPERTIES
%


% checks for an initialization function
if isfield(Sys, 'init_fun')
    P = Sys.init_fun(P);
end

if isfield(P, 'init_fun')
    P = P.init_fun(P0);
end

if (~isfield(Sys, 'type'))
    Sys.type = '';
end

if isfield(Sys, 'type')
    if strcmp(Sys.type,'traces')
        Pf = P;
        return ;
    end
end

if isfield(Sys,'init_fun')
    P = Sys.init_fun(P);
end
if isfield(P,'init_fun')
    P = P.init_fun(P);
end

if exist('is','var')
    org_dims = P.dim; % we save original dim and epsi
    org_epsi = P.epsi;
    P = SAddUncertainParam(P,is);
    P = SDelUncertainParam(P,org_dims,is);
end

if iscell(tspan)
    if(numel(tspan)==2)
        T = [tspan{1} tspan{2} tspan{2}];
    else
        T = cell2mat(tspan);
    end
else
    T = tspan;
end

InitSensi(Sys,P);

if isfield(P,'XS0') && isempty(P.XS0)
    P = rmfield(P,'XS0');
end

if ~isfield(P,'XS0')
    dims = P.dim;
    Ns = numel(dims);
    N = P.DimX;
    ix0 = dims(dims<=N); % Parameters in x0
    
    yS0 = zeros(N,Ns);
    
    for i=1:numel(ix0); % pour chaque condition initiale incertaine
        yS0(dims(i),i) = 1;
    end
    
    xS0 = reshape(yS0,N*Ns,1); % we stack columns of yS0
    
    P.XS0 = repmat(xS0,[1 size(P.pts,2)]);
    
end

Pf = cvm(93,P,T); %NM : what happens if P.XS0 is not empty?
CVodeFree();

if exist('is','var')
    Pf = SAddUncertainParam(Pf,org_dims);
    Pf = SDelUncertainParam(Pf,is,org_dims);
    Pf.epsi = org_epsi; %NM : not sure of that (wonder if cvm makes changes on epsi - I guess no)
end

end
