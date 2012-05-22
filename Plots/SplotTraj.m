function SplotTraj(S,proj,iX,opt,t0)%
%   SplotTraj(S [, proj, iX, opt, t0])
%
%   Plots the trajectories in field traj of S, recursively (explore children
%   structure, if any), in the phase space.
%   
%   Note:  Uses the plotting options defined in field traj_plot_opt, and
%   project on dimensions specified by field 'plot_proj' if these fields
%   are defined
%   
%   Prerequisite: S has a traj field. Else what's the point ?
%   
%   Inputs: 
%   
%    -  S        box sampling set. 
%    -  proj     variables to plot (all if [])
%    -  iX       indices of the initial pts from which to plot   
%    -  opt      plotting option e.g: {'r','LineWidth',4}
%    -  t0       starting time to plot trajectories from.
% 
%   Outputs:     
%      
%    -  Some figure, hopefully an interesting one.
%  
  
% Check inputs
  if (isfield(S, 'time_mult'))
    time_mult= S.time_mult;
  else  
    if (isfield(S,'opt'))
      if (isfield(S.opt, 'time_mult'))
        time_mult= S.opt.time_mult;
      else
        time_mult=1;
      end  
    else
      time_mult=1;
    end
  end   
  
   if (isfield(S, 'rescale'))
    rescale= S.rescale;
  else  
    if (isfield(S,'opt'))
      if (isfield(S.opt, 'rescale'))
        rescale= S.opt.rescale;
      else
        rescale=1;
      end  
    else
      rescale=1;
    end
  end   
  
  
  
    
  if (isempty(S.pts))
    error('S empty !');
    return
  end

  if (~isfield(S,'traj'))
    error('No trajectory computed for this set')
    return
  end

  if (~exist('t0'))
    t0= 0;
  end
   
  if (isfield(S,'traj_plot_opt'))
   
    opt = S.traj_plot_opt;
  elseif (~exist('opt')||isempty(opt))
    opt = {'b','LineWidth',1};
  end
  
  if (ischar(opt))
    opt = {opt};
  end   
      
  % find the projected axes
  
  if (isfield(S,'plot_proj'))
    proj = S.plot_proj;
  else 
    if (~exist('proj')||isempty(proj))       
      switch (S.DimX)
       case {1}
        proj=1;
       case {2}
        proj=[ 1 2 ];
       otherwise
        proj = [ 1 2 3 ];
      end
    end
  end

  if (~isnumeric(proj))
    stproj = proj;
    proj=[];
    for i = 1:numel(stproj)
      proj(i) = FindParam(S,stproj{i});
    end
  end  
  proj = proj(proj~=0);
  
  if (~exist('iX')||isempty(iX))
    iX = 1:size(S.pts,2);
  end
  
  switch (numel(proj))
    
   case {1}
    hold on;
    xlabel('time')
    
    if isfield(S,'ParamList')            
      ylabel(S.ParamList{proj(1)},'Interpreter','none');
    else
      ylabel(['x_' num2str(proj(1))]);
    end       
    
    nb = numel(S.traj);
           
    if (isfield(S,'traj_plot_opt'))
      opt = S.traj_plot_opt;
    end
      
    for i = iX
      time = S.traj(i).time;
      if (numel(time>t0))
        y = S.traj(i).X(proj(1),time>=t0) * rescale;
        x = S.traj(i).time(time>=t0)*time_mult;
        plot(x,y, opt{:});
        % drawnow
      end
    end
    drawnow       
      
   case {2}
    hold on;
    
    if isfield(S,'ParamList')            
      xlabel(S.ParamList{proj(1)},'Interpreter','none');
      ylabel(S.ParamList{proj(2)},'Interpreter','none');
    else
      xlabel(['x_' num2str(proj(1))]);
      ylabel(['x_' num2str(proj(2))]);
    end
       

    nb = numel(S.traj);
    
    if (isfield(S,'traj_plot_opt'))
      opt = S.traj_plot_opt;
    end
    
    for i = iX
      time = S.traj(i).time;
      if (numel(S.traj(i).time)>0)
        x = S.traj(i).X(proj(1),time>=t0)*rescale;
        y = S.traj(i).X(proj(2),time>=t0)*rescale;        
        plot(x,y,opt{:});        
      end 
    end     
    drawnow            
   
   otherwise
    hold on;
    
    if isfield(S,'ParamList')            
      xlabel(S.ParamList{proj(1)},'Interpreter','none');  
      ylabel(S.ParamList{proj(2)},'Interpreter','none');
      zlabel(S.ParamList{proj(3)},'Interpreter','none');
    else
      xlabel(['x_' num2str(proj(1))]);
      ylabel(['x_' num2str(proj(2))]);
      zlabel(['x_' num2str(proj(3))]);
    end
        
  
    nb = numel(S.traj);
    
    if (isfield(S,'traj_plot_opt'))
      opt = S.traj_plot_opt;
    end
    
    for i = iX
      time = S.traj(i).time;
      if (~isempty(S.traj(i).time))
        x = S.traj(i).X(proj(1),time>=t0)*rescale;
        y = S.traj(i).X(proj(2),time>=t0)*rescale;
        z = S.traj(i).X(proj(3),time>=t0)*rescale;
        plot3(x,y,z,opt{:});        
      end      
      drawnow
    end
  end
  
  grid on;
  if (t0==0)&&(numel(proj)>1)
    SplotPts(S,proj,iX);      
    if isfield(S,'Xf')      
      SplotXf(S,proj,iX);  
    end
  end
  