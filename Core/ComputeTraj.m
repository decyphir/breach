function Sf = ComputeTraj(Sys, S0,tspan)
%
%   Sf = ComputeTraj(Sys,S0,tspan) 
%  
%   Compute trajectories issued from points in S0 on the
%   time interval tspan
%    
%   Inputs: 
%    
%    -  Sys      System (needs to be compiled)
%    -  S0       Initial sampling of the uncertain set  
%    -  tspan    interval of the form [t0, tf]; Fixed time instants can
%    also be specified tspan = [t0 t1 ... tN];
%
%   Outputs:
%      
%    -  Sf       Sampling structure augmented with the field traj
%                containing computed trajectories
%
 
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
    