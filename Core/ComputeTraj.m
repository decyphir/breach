function Sf = ComputeTraj(Sys, S0,tspan, u)
%  COMPUTETRAJ compute trajectories for a system given initial conditions
%  and parameters
%    
%  Synopsis:   Pf = ComputeTraj(Sys,P0,tspan [,u]) 
%  
%   Compute trajectories issued from points in S0 on the
%   time interval tspan
%    
%   Inputs: 
%    
%    -  Sys      System (needs to be compiled)
%    -  P0       Initial conditions and params given in a parameter set 
%                or in an array of size DimX x nb_traj
%    -  tspan    interval of the form [t0, tf]; Fixed time instants can
%                also be specified tspan = [t0 t1 ... tN];
%    -  u        is a structure array of nb_traj inputs. u(i) must be a
%                struct with fields 
%                    - params_idx : indicates which parameters are made 
%                        time dependant  
%                    - tin : times when the input changes    
%                    - values : values of the parameters       
%     
%   Outputs:
%      
%    -  Sf       Sampling structure augmented with the field traj containing
%                 computed trajectories if the input is a param set 
%                
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
 
  if (isfield(Sys, 'type'))

    % If system type is only traces, check consistency of params and pts
    if strcmp(Sys.type, 'traces')
      Sf=S0;
      for (i = 1:numel(Sf.traj))
        Sf.traj(i).param = Sf.pts(:,i)';      
      end
      return;
    end
  end

  if (isfield(S0,'traj'))
    if ((numel(S0.traj(1).time)==numel(tspan)))
      if (sum(S0.traj(1).time==tspan))
        Sf = S0;
        return;
      end
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
  
  
  if (exist('u'))
     
    err = check_u(u);
    if (err~=0)
      error(err);
    end
    
    Sf=cvm(61,S0,T,u);            
    Sf.u = u;
    
  else
    Sf=cvm(61,S0,T);        
  end    
  
  CVodeFree();
  
  if output_trajs
    Sf = Sf.traj;
  end
  
function err = check_u(u)
  
  err = 0;

  if ~isstruct(u)
    err = 'u has to be a structure';
    return;
  end
    
  if ~isfield(u,'params_idx')
    err = 'missing field params_idx';
    return;
  end
  
  if ~isfield(u,'time')
    err = 'missing field time';
    return;
  end
  
  if ~isfield(u,'values')
    err = 'missing field values';
    return;
  end
    
  if numel(u.params_idx) ~= size(u.values, 1)
    err = 'numel(u.params_idx) should be equal to  size(u.values, 1)';
    return;
  end
  
  if  numel(u.time) ~= size(u.values, 2)
    err = 'numel(u.time) should be equal to size(u.values, 2)'; 
    return;
  end
    
