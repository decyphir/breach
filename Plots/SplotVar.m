function SplotVar(S,iX,ipts,opt, bool_same_axe)
% SPLOTVAR Plots trajectories variables separatly
%
% Synopsis:  SplotVar(P, [iX,ipts,opt])
%     
%   Inputs: 
%   
%    -  P        Parameter set 
%    -  iX       indices of the  X variables to plot (optional, absent
%    means all)
%    -  ipts     indices of the trajectories in S to plot (optional,
%    absent means all)
%    -  opt      plotting options    
%    - bool_same_axe  if 1, all variables are plot on the same axe
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
    return;
  end
  
  
  if (~exist('iX')||isempty(iX))
    iX = 1:S.DimX;
  end
  
  if (~exist('ipts')||isempty(ipts))
    ipts = 1:numel(S.traj);
  end
  
  if ( isfield(S, 'traj_ref') )    
    ipts = unique(S.traj_ref(ipts));        
  end
     
  if (~isempty(iX))
    if (~isnumeric(iX))  
      if isstr(iX)
          iX = {iX};
      end
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
  

  if (~exist('bool_same_axe'))
    same_axe = 0;
  else
    same_axe = bool_same_axe;
  end
  
  if (isfield(S,'plot_proj'))
    proj = S.plot_proj;
  else
    proj = 1:size(S.pts,1);
  end
     
  
  if (nargin == 2)
    ipts = 1:numel(S.traj);
  end
  
  
  % Plot options
  
  colors = hsv(numel(ipts));
  colors = colors(:,[3 2 1]);
  
  if (~exist('opt')||isempty(opt))
    if (isfield(S,'traj_plot_opt'))
      opt = S.traj_plot_opt;
    else
      opt = [];
   end
  end
     
  
  if (same_axe==1)
    lg = {};
    for i = ipts
      
      time = S.traj(i).time;       
      grid on;
      % set(gca,'FontSize',12,'FontName','times');
      hold on;  
                  
      X = S.traj(i).X(iX(:),:);       
      plot(time*time_mult,X);
            
    end
    
    for j = 1:numel(iX)
      if isfield(S,'ParamList')            
        lg = {lg{:} S.ParamList{iX(j)} }; 
        
      else
        lg = {lg{:} ['x_' num2str(iX(j))] }; 
        
      end
    end
    hl = legend(lg);
    set(hl, 'Interpreter','none');
    hold off;
    xlabel('time')
  
  else % plots on multi axes       
    ci=1;
    for i = ipts
      
      time = S.traj(i).time;       
      
      for j = 1:numel(iX)
        if (numel(iX)>1)
          subplot(numel(iX),1,j)
        end
        grid on;
        %    set(gca,'FontSize',12,'FontName','times');
        hold on;  
        
        if isfield(S,'ParamList')            
          ylabel(S.ParamList{iX(j)},'Interpreter','none');
        else
          ylabel(['x_' num2str(iX(j))]);
        end
        
        x = S.traj(i).X(iX(j),:);       
        if isempty(opt)
            plot(time*time_mult,x,'Color', colors(ci,:));   
        else
          plot(time*time_mult,x,opt{:});
        end
      end
    ci = ci+1;  
    end
    hold off;
    xlabel('time')  
  end

  
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

