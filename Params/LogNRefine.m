function P = LogNRefine(P, N)
%LOGNREFINE refines randomly using log normal distribution. epsi is
% interpreted as standard deviation
%
% Synopsis: P = LOGNREFINE(P, N)
%
% Input:
%  - P : the parameter set to refine
%  - N : The number of parameter set to generate for each set of parameter
%        values in P.
%
% Output:
%  - P : The new parameter set
%
%TODO : THIS FUNCTION DOES NOT MANAGE THE TRAJ_REF AND TRAJ_TO_COMPUTE
% FIELDS
%
% Example (Lorentz84):
%
% TO WRITE
%
%See also QuasiRefine RandomLogRefine Refine
%

if(N<=1)
    return;
end

dim_num = numel(P.dim);
pts = kron(P.pts,ones(1,N));

for ii = 1:dim_num
    for jj = 1:size(P.pts,2)
        m = P.pts(P.dim(ii),jj)^2;
        v = P.epsi(ii,jj)^2;
        mu = log(m/sqrt(v+m));
        sigma = sqrt(log(v/m+1));
        pts(P.dim(ii),1+N*(jj-1):N*jj) = lognrnd(mu,sigma,1,N); %NM : improve by using the matricial version
    end
end
P.pts = pts;
%P.epsi = kron(P.epsi,ones(1,size(P.pts,2)))/(floor(N^(1/dim_num)));
P.epsi = kron(P.epsi,ones(1,size(P.pts,2)))/(N^(1/dim_num));

if isfield(P,'selected')
    P.selected = zeros(1, size(P.pts,2));
end

P = Preset_traj_ref(P);

end
