function Sr = LOGNREFINE(S, N)
% LOGNREFINE  refine randomly using log normal distribution  
%
% usage: Sr = LOGNREFINE(S0, N)
%
% Interprets epsi as standard deviation
%  

  dim_num = numel(S.dim);
  Sr = S;
  Sr.pts = repmat(S.pts,[1 N]);
  Sr.epsi = kron(Sr.epsi,ones(1,size(Sr.pts,2)))/(floor(N^(1/dim_num)));

  for i = 1:numel(S.dim)
  
    m = S.pts(S.dim(i),1);
    v = S.epsi(i,1)^2;
    mu = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1)); 
    Sr.pts(S.dim(i),:)  = lognrnd(mu,sigma,1,N);
    
  end
  
 
  if (isfield(S,'selected'))
    Sr.selected = zeros(1, size(S.pts,2));    
  end