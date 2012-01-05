function Sf = ComputeTraj(Sys, S0,tspan)
%  COMPUTETRAJ compute trajectories for a system given initial conditions
%  and parameters
%    
%  Synopsis: 
%  Sf = ComputeTraj(Sys,S0,tspan) or  trajs = ComputeTraj(Sys, P0,tspan)  
%  
%   Compute trajectories issued from points in S0 on the
%   time interval tspan
%    
%   Inputs: 
%    
%    -  Sys      System (needs to be compiled)
%    -  S0       Initial conditions and params given in a parameter set 
%                or in an array of size DimX x nb_traj
%    -  tspan    interval of the form [t0, tf]; Fixed time instants can
%    also be specified tspan = [t0 t1 ... tN];
%
%   Outputs:
%      
%    -  Sf       Sampling structure augmented with the field traj
%                containing computed trajectories if the input is a param
%                set 
% or - trajs     array of trajectories
%  
 
  
% checks if we have a parameter set or 
    
  output_trajs = 0;  
  if ~isstruct(S0)
    
    if (size(S0,1) ~= Sys.DimP)
      if (size(S0,2) == Sys.DimP) % be smart, try transpose in case it works
        S0 = S0';
      else
        error('Second argument must be a parameter set or be of dimension Sys.DimP x ?')
      end
    end
    output_trajs = 1;
    
    pts = S0;
    S0 = CreateSampling(Sys,1);
    S0.pts = pts;
    S0.epsi = ones(1, size(pts,2));

  end
 
  
  if (isfield(S0,'traj'))
    if ((numel(S0.traj(1).time)==numel(tspan)))
      if (sum(S0.traj(1).time==tspan))
        Sf = S0;
        return;
      end
    end
  end
  
  if (isfield(Sys, 'type'))
    if strcmp(Sys.type, 'traces')
      Sf=S0;
      return;
    end
  end
 
  InitSystem(Sys);

  if iscell(tspan)
    if (numel(tspan)==2)
      T = [tspan{1} tspan{2} tspan{2}];  % not really nice.. should be
                                         % changed some day
    else
      T= cell2mat(tspan);
    end
  else
    T = tspan;
  end
  
  if (isfield(S0,'test'))
    switch S0.test
     case 1
      Sf=cvm(61,S0,T,1);            
    end
  else
    Sf=cvm(61,S0,T);        
  end    
  
  CVodeFree();
  
  if output_trajs
    Sf = Sf.traj;
  end
  