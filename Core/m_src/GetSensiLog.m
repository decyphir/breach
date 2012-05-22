function dlx = GetSensiLog(traj,i,is) 
%
%   dlxx = GetSensiLog(traj,i, is)
%  
%   Returns the logarithmic sensitivity of xi with respect to parameter
%   is
%
  
  % get  dX/dp[t]        
  
  dx = traj(i).XS(is,:);    

  % get X[t]        
  x = traj(i).X(iX(j),:); 
  
  % replace zeros by small quantities
  ind = find(abs(x)<1e-16);  
  x(ind) = sign(x(ind))*1e-16;       
  x(x==0) = 1e-16;

  % get p
  p = traj(i).param(is);    
        
  dlx=  (dx*p)./x;
