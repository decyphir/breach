function SplotVarExp(S,iX,ipts)
%
%   SplotVarExp(S, [iX,ipts])
%
%   Plots trajectories states variables and represent expansion given by
%   sensitivity
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
%
%   Outputs:     
%      
%    -  Some figure, hopefully an interesting one.
%  
  
  % Check inputs
  if (isempty(S.pts))
    disp('S empty !');
    return
  end

  if (~isfield(S,'traj'))
    disp('No trajectory computed for this set')
  end
    
  if (nargin == 1)
    iX = 1:S.DimX;
    ipts = 1:numel(S.traj);
  elseif (nargin == 2)
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
  
  
  depth=1;
  Sp = S;
  color_rez=0;
  
  while(isfield(Sp,'child'))
    Sp = Sp.child;
    depth = depth+1;
  end
  
  if (isfield(S,'traj_plot_opt'))
    color_rez = 0;
    opt = S.traj_plot_opt;
  else
    opt = {'b','LineWidth',1};
  end

  colors = gray(depth);
  colors(:,[1 2]) =0;
  colors = colors(end:-1:1,:);

  if (depth==1)
    colors = [0 0 1];
  end   
  
  Sp=S;

  if (isfield(S,'plot_proj'))
    proj = S.plot_proj;
  else
    proj = 1:size(S.pts,1);
  end
     
    for k=1:depth
      
      if (isfield(Sp,'traj_plot_opt'))
        opt = Sp.traj_plot_opt;
      end
      

      if (nargin == 2)
        ipts = 1:numel(Sp.traj);
      end
      
      for i = ipts
        
        time = Sp.traj(i).time;       
        
        for j = 1:numel(iX)
          subplot(numel(iX),1,j)
          hold on;  
         
          if isfield(S,'ParamList')            
            ylabel(S.ParamList{iX(j)},'Interpreter','none');
          else
            ylabel(['x_' num2str(iX(j))]);
          end
         
          x = Sp.traj(i).X(iX(j),:);       
          e = Sp.traj(i).Expa(iX(j),:);          
          base = min(x-e-.1);
          
          if (color_rez)
            plot(time, x,opt{:},'Color',colors(k,:));
          else
            
            area(time,x+e,base,'FaceColor',[.5 .9 .6],...
                 'EdgeColor','w',...
                 'LineWidth',1)
            hold on
            area(time,x-e,base,'FaceColor',[1 1 1],...
                 'EdgeColor','w',...
                 'LineWidth',1)
            plot(time,x+e,'k');
            plot(time,x-e,'k');
            plot(time,x,opt{:});
          end        
        end
      end      
      xlabel('time')
      if (k<depth)
        Sp = Sp.child;
      end
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

  