function Sf = ComputeTrajSensi(Sys,S,tspan)
%
%   Sf = ComputeTrajSensi(Sys,S0,tspan) 
%  
%   Compute trajectories with corresponding sensitivities issued from points
%   in S0 on the time interval tspan.
%    
%   Inputs: 
%   
%    -  Sys      System (needs to be compiled)  
%    -  S0       Initial sampling of the uncertain set  
%    -  tspan    interval of the form [t0, tf];
%
%   Outputs:
%      
%    -  Sf       Sampling structure augmented with the field traj
%                containing computed trajectories
%
  

  if (isfield(Sys, 'fake'))
    Sf=S0;
    return;
  end  

  if iscell(tspan)
      if (numel(tspan)==2)
        T = [tspan{1} tspan{2} tspan{2}];
      else
        T= cell2mat(tspan);
      end
    else
      T = tspan;
    end
     
  InitSensi(Sys,S);
  
  if (isfield(S,'XS0'))
    if (isempty(S.XS0))
      S = rmfield(S,'XS0');
    end
  end
      
  if (~isfield(S,'XS0'))
    dims = S.dim;
    Ns = numel(dims);
    N = S.DimX;
    ix0 = dims(dims<=N); % Parameters in x0
    ip = dims(dims>N); % Parameters NOT in x0
    
    xS0=[];
    yS0 = zeros(N,Ns);
    
    for i=1:numel(ix0);
      yS0(dims(i),i)= 1;
    end
    
    for i=1:Ns
      xS0 = [xS0 ; yS0(:,i)];
    end
    
    S.XS0 = repmat(xS0,[1 size(S.pts,2)]);
      
  end

  Sf = cvm(93,S,T);
  CVodeFree();
    
  