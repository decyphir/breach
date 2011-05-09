function SplotVar(S,iX,ipts,opt)
%
%   SplotVar(S, [iX,ipts,opt])
%
%   Plots trajectories variables separatly
%   
%   Note:   Uses the plotting options defined in field traj_plot_opt, and
%   project on dimensions specified by field 'plot_proj'.
%   
%   Prerequisite: S has a traj field. Else what's the point ?
%   
%   Inputs: 
%   
%    -  S        rectangular sampling set. 
%    -  iX       indices of the  X variables to plot (optional, absent
%    means all)
%    -  ipts     indices of the trajectories in S to plot (optional,
%    absent means all)
%    -  opt      plotting options    
%  
%   Outputs:     
%      
%    -  Some figure, hopefully an interesting one.
%  

  % Check inputs


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
    disp('No trajectory computed for this set')
  end

  if (~exist('iX')||isempty(iX))
    iX = 1:S.DimX;
  end
  
  if (~exist('ipts')||isempty(ipts))
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
  
  if (~exist('opt'))
    if (isfield(S,'traj_plot_opt'))
      opt = S.traj_plot_opt;
    else
      opt = {'b','LineWidth',1};
    end
  end

  if (isfield(S,'plot_proj'))
    proj = S.plot_proj;
  else
    proj = 1:size(S.pts,1);
  end
     
  
  if (nargin == 2)
    ipts = 1:numel(S.traj);
  end
  
  for i = ipts
    
    time = S.traj(i).time;       
    
    for j = 1:numel(iX)
      if (numel(iX)>1)
        subplot(numel(iX),1,j)
      end
      grid on;
      set(gca,'FontSize',12,'FontName','times');
      hold on;  
      
      if isfield(S,'ParamList')            
        ylabel(S.ParamList{iX(j)});
      else
        ylabel(['x_' num2str(iX(j))]);
      end
      
      x = S.traj(i).X(iX(j),:);       
      
      plot(time*time_mult,x,opt{:});
    end
  end      
  hold off;
  xlabel('time')
  

  
function index=  FindParam(Sys,param)
  
  if ~isfield(Sys,'ParamList')
    error('No parameter list ...');
  end

  for j = 1:numel(Sys.ParamList)
    if (strcmp(Sys.ParamList{j}, param));
      index = j;
      return;            
    end
  end
  
  error(['Parameter ' param ' not found']);

