function SplotETraj(S,t0,d)

%   SplotETraj(S [,t0])
%
%   Plots the trajectories in field etraj, recursively (explore children structure).
%   
%   Note:   Uses the plotting options defined in field traj_plot_opt, and
%   project on dimensions specified by field 'plot_proj'.
%   
%   Prerequisite: S has an etraj field. Else what's the point ?
%   
%   Inputs: 
%   
%    -  S        rectangular sampling set. Any set with a field 'Pts'
%                will do though.
%    -  t0       starting time to plot traj from.
% 
%   Outputs:     
%      
%    -  Some figure, hopefully an interesting one.
%  
  
  % Check inputs
  if (isempty(S.pts))
    error('S empty !');
    return
  end

  if (~isfield(S,'etraj'))
    error('No estimated trajectory computed for this set')
    return
  end

  if (nargin==1)
    t0 = 0;
    d =0;
  elseif (nargin == 2)
    d=0;
  end

  depth=1;
  Sp = S;
  color_rez=1;
  
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
  
  if (~exist('iX')||isempty(iX))
    iX = 1:size(S.pts,2);
  end
  
  
  switch (numel(proj))
  
   case {1}
    hold on;
    xlabel('time')
    ylabel(['x_' num2str(proj(1))]);
       
    for k=1:depth
      nb = numel(Sp.etraj);
           
      if (isfield(Sp,'traj_plot_opt'))
        opt = Sp.traj_plot_opt;
      end
      
      for i = 1:nb
        time = Sp.etraj(i).time;
        if (numel(time>t0))
          x = Sp.etraj(i).X(proj(1),time>=t0);
          y = Sp.etraj(i).time(time>=t0);
          plot(x, y,opt{:})          
          % drawnow
        end
      end
      
      if (k<depth)
        Sp = Sp.child;
      end      
      drawnow       
    end
   
   case {2}
    hold on;

    xlabel(['x_' num2str(proj(1))]);
    ylabel(['x_' num2str(proj(2))]);

    for k=1:d
      if isfield('child',Sp)
        Sp = Sp.child;
      end
    end
    
    for k=d+1:depth 
       nb = numel(Sp.etraj);

       if (isfield(Sp,'traj_plot_opt'))
         opt = Sp.traj_plot_opt;
       end
       for i = 1:nb
         time = Sp.etraj(i).time;
         if (numel(Sp.etraj(i).time)>0)
           x = Sp.etraj(i).X(proj(1),time>=t0);
           y = Sp.etraj(i).X(proj(2),time>=t0);
           if (color_rez)
             plot(x, y,opt{:},'Color',colors(k,:));
           else
             plot(x,y,opt{:});
           end
         end 
       end
     
       if (k<depth)
         Sp = Sp.child;
       end
       drawnow       
     end
       
   otherwise
    hold on;
    
    xlabel(['x_' num2str(proj(1))]);
    ylabel(['x_' num2str(proj(2))]);
    zlabel(['x_' num2str(proj(3))]);
    
    for k=1:d
      if isfield('child',Sp)
        Sp = Sp.child;
      end
    end
    
    for k=d+1:depth
      nb = numel(Sp.etraj);

      if (isfield(Sp,'traj_plot_opt'))
        opt = Sp.traj_plot_opt;
      end

      for i = 1:nb
        time = Sp.etraj(i).time;
        if (~isempty(Sp.etraj(i).time))
          x = Sp.etraj(i).X(proj(1),time>=t0);
          y = Sp.etraj(i).X(proj(2),time>=t0);
          z = Sp.etraj(i).X(proj(3),time>=t0);
       
          if (color_rez)
            plot3(x, y,z,opt{:},'Color',colors(k,:));
          else
            plot3(x,y,z,opt{:});
          end
        end
      end
       
      drawnow
            
      if (k<depth)
        Sp = Sp.child;
      end
    end
  end

  grid on;
  if (t0==0)
    SplotPts(S);      
  end
  
