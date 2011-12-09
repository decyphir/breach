function SplotSensi(S,iX,iP,ipts)
%
%   SplotSensi(S,iX,iP,ipt)
%
%   Plots trajectories sensibilities of state variables iX
%   w.r.t. parameters iP
%   
%   Note: Uses the plotting options defined in field traj_plot_opt, and
%   project on dimensions specified by field 'plot_proj'.
%   
%   Prerequisite: S has a traj field with an XS field. Else what's the point ?
%   
%   Inputs: 
%   
%    -  S        rectangular sampling set. 
%
%    -  iX       indice of the X variables for which to plot the sensitivity (idem)
%   
%    -  iP       indice of the sensitive parameter (optional, all if absent)
%
%    -  ipts     (optional) indices of the trajectories in S to plot (all
%                if absent)
%
%   Outputs:     
%      
%    -  Some figure, hopefully an interesting one.
%  
  
%   Check inputs 
     
  
  if (isfield(S, 'time_mult'))
    time_mult= S.time_mult;
  else
    time_mult=1;
  end
  
  
  if (isempty(S.pts))
    disp('S empty !');
    return
  end
  
  if (~isfield(S,'traj'))
    disp('No trajectory computed for this sampling')
  end
  
  if (nargin == 1)
    iX = 1:S.DimX;
    iP =S.dim;
    ipts = 1:numel(S.traj);
  elseif (nargin == 2)
    iP = S.dim;
    ipts = 1:numel(S.traj);
  elseif (nargin == 3)
    ipts = 1:numel(S.traj);
  end

 if (~isempty(iX))
    if (~isnumeric(iX))  
      NiX = iX;
      iX = [];
      for i = 1:numel(NiX)
        ind = FindParam(S,NiX{i});
        iX(i) = ind;
      end    
    end
  else
    iX = 1:S.DimX;
  end
  iX = iX(iX<=S.DimX);

  if (isempty(iP))
    iP = S.dim;    
  elseif (~isnumeric(iP))  
    NiP = iP;
    iP = [];
    for i = 1:numel(NiP)
      ind = FindParam(S,NiP{i});
      iP(i) = ind;
    end    
  end
    
 
  if (isfield(S,'traj_plot_opt'))
    opt = S.traj_plot_opt;
  else
    opt = {'b','LineWidth',1};
  end

  colors = hsv(numel(iP));

  if (isfield(S,'plot_proj'))
    proj = S.plot_proj;
  else
    proj = 1:size(S.pts,1);
  end
       
  if (isfield(S,'traj_plot_opt'))
    opt = S.traj_plot_opt;
  end
          
  for i = ipts
    
    time = S.traj(i).time;       

    for j = 1:numel(iX)
      if (numel(iX)>1)
        subplot(numel(iX),1,j)
      end
      hold on;
      ylabel(['Sensi(' S.ParamList{iX(j)} ')'],'Interpreter','none');

      leg = {};

      for k = 1:numel(iP)
        is = (find(S.dim==iP(k))-1)*size(S.traj(i).X,1)+iX(j);      
        x = S.traj(i).XS(is,:);  
        plot(time*time_mult, x,opt{:},'Color',colors(k,:));
        leg = {leg{:}, S.ParamList{iP(k)}};
      
      end
      legend(leg{:});
    end
  end
  
  xlabel('time','Interpreter','none');
  
  end
  
  
  