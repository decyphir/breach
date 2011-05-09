function [Sf ErrStat] = sreach(Sys,S0,tspan,options)
%
%   [Sf Err] = sreach(Sys,S0,tspan,options)
%
%   Reachability routine using sensitivity.   
%    
%   Inputs: 
%
%    -  Sys      System  
%    -  S0       Initial parameter set  
%    -  tspan    Timespan
%  
%    - options    can contain the following fields:   
%
%       .tol      absolute tolerance vector (or scalar) for nonlin error
%       .delta     vector of dimension n, where n is the dimension of the
%                  system; The algorithm won't refine S0 into boxes of
%                  dimensions less than delta  
%       .NbIterErr maximum number of refinement iterations  
%  
%   Outputs:
%      
%    -  Sf       Sampling of the reachable set
%    -  Err      Error estimation history
%   

try
   S0= rmfield(S0, 'selected');
end
  dimx = Sys.DimX;
  
  if (nargin == 3)
    options = [];
  end
  
  if (isfield(options,'tol'))
    tol = options.tol;
    if (size(tol,1)==1)
      tol= tol';
    end
  else
    tol = .1;
  end
    
  if (isfield(options,'NbIterErr'))
    NbIterError = options.NbIterErr;
  else
    NbIterError = 5;
  end  
      
  if (isscalar(tol))
    tol = repmat(tol,dimx,1);
  end  

  if isfield(options,'delta')
    if ~isempty(options.delta)
      delta = options.delta;
    else
      min_epsi = min(S0.epsi);
      [ m ind ]=  min(min_epsi);
      delta = S0.epsi(:,ind)/2^NbIterError;
    end
    
  else
    min_epsi = min(S0.epsi);
    [ m ind ]=  min(min_epsi);
    delta = S0.epsi(:,ind)/2^NbIterError;
  end
  
  
  if (size(delta,1) ==1)
    delta = delta';
  end
  
  % Refine sampling until sensitivity provides a good approximation of the
  % reachable set
  
  nbiter = 0;
  global nbtraj_computed;
  nbtraj_computed = 0;
  SToCheck = ComputeTrajSensi(Sys,S0, tspan);
  nbtraj_computed = nbtraj_computed+numel(SToCheck.traj);
  Sf = SCreateNull(S0);
  Sf = rmfield(Sf,'traj');
  Snotr = Sf;
  ErrStat = []; 
  
  while (numel(SToCheck.traj))
    
    Dp = repmat(delta,[1 size(SToCheck.epsi,2)]);
    kn = find(sum(SToCheck.epsi>Dp,1));
    nkn = find(prod(1*(SToCheck.epsi<=Dp),1));
    
    Snotrn = select(SToCheck,nkn);    
    SToCheck = select(SToCheck,kn);
    
    if (isfield(Snotrn,'ToCheck'))
      Snotrn = rmfield(Snotrn,'ToCheck');
      if (size(Snotrn.pts,2))
        Snotr= SConcat(Snotr,Snotrn);
      end
    end
    if (numel(kn)==0)
      break;
    end
   
    [Sfn SToCheck Err] = SCheckSensiEstim(Sys,SToCheck,tspan, tol);

    ErrStat = [ErrStat; Err];
    Affiche(Err);
    if (size(Sfn.pts,2))
      Sf = SConcat(Sf,Sfn);
    end
    nbiter = nbiter+1;
    if nbiter>=NbIterError
      break
    end
  end
  
  Sf= SConcat(Sf,Snotr);
  Sf.selected = zeros(1,size(Sf.pts,2));
  
  if (size(SToCheck.pts,2))
    if (isfield(SToCheck,'ToCheck'))
        SToCheck.selected = SToCheck.ToCheck;
        SToCheck = rmfield(SToCheck,'ToCheck');
    else
        SToCheck.selected = zeros(1,size(SToCheck.pts,2));
    end
    Sf = SConcat(Sf,SToCheck);
  end
 
%  nbtraj_computed
  
% $$$ function S = Select(S, kn)              
% $$$   
% $$$   S.pts = S.pts(:,kn);
% $$$   S.epsi = S.epsi(:,kn);
% $$$   S.Xf =  S.Xf(:,kn);
% $$$   S.ExpaMax = S.ExpaMax(:,kn);
% $$$   S.XSf = S.XSf(:,kn);        
% $$$  
%  S.traj = S.traj(:,kn);
  
function Affiche(Err)
  
  fprintf('Error:');
  for i = 1:numel(Err)
    fprintf(' %g', Err(i));
  end
  fprintf('  ');
