function Y = GetTrajValues(Pf, iX,  t)
% GETTRAJVALUES  extract values of a variable evolution in a set of trajectories 
%
% Synopsis: Y = GetTrajValues(Pf,iX,t)
%    
%  Pf : parameter set with pre-computed trajectories
%  iX : index or name of the required variable
%  t  : time instant when to extract values
%
  
  if (~isfield(Pf,'traj'))
    error('Compute trajectories first')
  end
  
  if (ischar(iX))
    iX = FindParam(P, iX);
  end
  
  
  X = cat(1, Pf.traj.X);
  
  X = X(iX:Pf.DimX:end,:);   
  
  Y = interp1(Pf.traj(1).time, X',t)';