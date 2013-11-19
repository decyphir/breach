function Pf = ComputeTrajExp(Sys, P, tspan)
%COMPUTETRAJEXP computes maximum expansion of trajectories issued from
% points in (the root of) P on the time interval tspan. Does not return the
% trajectories themselves.
%
% Synopsis: Sf = ComputeTrajExp(Sys,P,tspan)
%
% Inputs:
%  -  Sys   : System (needs to be compiled)
%  -  P     : Initial parameter set
%  -  tspan : interval of the form [t0, tf];
%
% Outputs:
%  -  Pf : Sampling structure augmented with the field traj containing
%          computed trajectories
%

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

if(isfield(P,'XS0') && isempty(P.XS0))
    P = rmfield(P,'XS0');
end

if ~isfield(P,'XS0')
    
    if(tspan(1) == 0.)
        dims = P.dim;
        Ns = numel(dims);
        N = P.DimX;
        ix0 = dims(dims<=N); % Parameters in x0
        ip = dims(dims>N); % Parameters NOT in x0
        
        nb_traj= size(P.pts,2);
        xS0 = [];
        yS0 = zeros(N,Ns);
        nb_traj = size(P.pts,2);
        
        for ii=1:numel(ix0);
            yS0(dims(ii),ii) = 1;
        end
        
        for ii=1:Ns
            xS0 = [xS0 ; yS0(:,ii)];
        end
        P.XS0 = repmat(xS0,[1 nb_traj]);
        
        
        if(tspan(2) == 0.)
            Pf = P;
            Pf.Xf = P.pts;
            Pf.XSf = P.XS0;
            Pf.tf = tspan(2)*ones(1,nb_traj);
            
            Expa0 = zeros(size(P.pts,1),1);
            Pf.traj = [];
            
            for ii = 1:nb_traj
                traj.time = tspan;
                traj.X = Pf.pts(:,ii);
                Expa0(P.dim) = Pf.epsi(:,ii);
                traj.Expa = [ Expa0  abs(Pf.XSf)*Pf.epsi(:,ii)];
                traj.U = [];
                Pf.traj = [Pf.traj traj];
            end
            
            return
            
        end
    end
end

Pf = cvm(92,P,tspan);

end
