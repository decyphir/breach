function Sf = ComputeExpansion(Sys,S,tspan)

%
%   Sf = ComputeExpansion(Sys,S0,tspan) 
%  
%   Compute maximum expansion of trajectories issued from points in (the root of) S0 on the
%   time interval tspan. Does not return the trajectories themselves.
%    
%   Prerequisite: 
%  
%       Some simulator must be correctly initialized, and
%       sensibilities must be activated.
%
%   Inputs: 
%   
%    -  S0       Initial sampling of the uncertain set  
%    -  tspan    interval of the form [t0, tf];
%
%   Outputs:
%      
%    -  Sf       Sampling structure augmented with the field traj
%                containing computed trajectories
%

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
    if (tspan(1) == 0.)
      
      dims = S.dim;
      Ns = numel(dims);
      N = S.DimX;
      ix0 = dims(dims<=N); % Parameters in x0
      ip = dims(dims>N); % Parameters NOT in x0
             
      nb_traj= size(S.pts,2);
      xS0=[];
      yS0 = zeros(N,Ns);
      
      for i=1:numel(ix0);
        yS0(dims(i),i)= 1;
      end
     

      for i=1:Ns
        xS0 = [xS0 ; yS0(:,i)];
      end
      S.XS0 = repmat(xS0,[1 nb_traj]);       
     
      if (tspan(2) == 0.)
        Sf= S;
        Sf.Xf = S.pts;
        Sf.XSf= S.XS0;
        Sf.tf = tspan(2)*ones(1,nb_traj);
    
        Expa0 = zeros(dimx,1);
        Sf.traj = [];
      
        for i = 1:nb_traj
          traj.time = tspan;
          traj.X = Sf.pts(:,i);
          Expa0(S.dim) = Sf.epsi(:,i);
          traj.Expa = [ Expa0  abs(Sf.XSf)*Sf.epsi(:,i)]; 
          traj.U = [];
          Sf.traj =[Sf.traj traj];
        end

        return 
      end
    end
  end
  
  
  Sf = cvm(90,S,tspan);
  
