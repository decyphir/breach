function Sf = ComputeTrajSensi(Sys,S,tspan, is)
%COMPUTETRAJSENSI compute trajectories with  sensitivities
%  
%  Synopsis:  Pf = ComputeTrajSensi(Sys,P0,tspan [, is]) 
%  
%   Compute trajectories with corresponding sensitivities issued from points
%   in P0 on the time interval tspan.
%    
%   Inputs: 
%   
%    -  Sys      System (needs to be compiled)  
%    -  P0       Initial parameter set
%    -  tspan    interval of the form [t0, tf], t0:dt:tf, etc
%    -  is       Parameter sensitivities to compute, if absent uses uncertain parameters in P
%  
%   Outputs:
%      
%    -  Pf       Sampling structure augmented with the field traj
%                containing computed trajectories with sensitivities
%
  
  if (exist('is'))  
    org_dims = S.dim;
    S = SAddUncertainParam(S,is);        
  end
    
  if (isfield(Sys, 'type'))
    if (Sys.type == 'traces')
      Sf=S;
      return;
    end
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
     
  if (exist('is'))    
    Sf = SDelUncertainParam(Sf,is,org_dims);  
  end