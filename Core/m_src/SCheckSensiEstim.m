function  [S SToCheckNew ErrMax] = SCheckSensiEstim(Sys,SToCheck, tspan,tol)
%
%  [S SToCheckNew ErrMax] = SCheckSensiEstim(Sys,SToCheck, tspan, tol)
%  
%  Check for each points x in S0.pts whether the estimation using
%  sensitivity matrix is good enough. If it is not, it replaces x by
%  refining hierarchically around x and marks that the new points need to be
%  checked. Otherwise, marks x as ok.
%  
  global nbtraj_computed;
  ErrMax = zeros(1,Sys.DimX);
  
  if (~isfield(SToCheck,'ToCheck'))
    SToCheck.ToCheck = ones(1, numel(SToCheck.traj));
  end

  St = RefineEstim(SToCheck,2);
  
  fprintf('still %d  traj to check ...\n', size(St.pts,2));  

  St = ComputeTrajSensi(Sys,St, tspan);    
  nbtraj_computed = nbtraj_computed+numel(St.traj);
  St = TrajErr(St);  
  St.ToCheck =[];
            
  for j = 1:numel(St.traj)
    
    l = numel(St.traj{j}.time);
    TOL = repmat(tol',l,1);        
    ExpaErrMax = max((St.traj{j}.ExpaErr')./TOL);

    %  ExpaErrMax= sqrt(sum((St.traj{j}.ExpaErr')./(ATOL+RTOL.*St.traj{j}.Expa').^2,2));
    ErrMax = max([ErrMax ; ExpaErrMax]);
        
    Err = max(ExpaErrMax);
    St.ToCheck(j)= (Err>1); % Check if the refined point is ok        
  end
 
  kToCheck = find(St.ToCheck);
  kNotToCheck = find(~St.ToCheck);
  
  S = Sselect(St,kNotToCheck);
  S = rmfield(S,'traj');
  S = rmfield(S,'etraj');
  SToCheckNew = Sselect(St, kToCheck);
