function SplotBoxTraj(S,proj,dt,c,alph, t0)
%
%  SplotBoxTraj(S,proj,tspan,c,alpha,t0)
% 
%
    
  if (~exist('c')||isempty(c))
    c = 'r';
  end
  if (~exist('dt')||isempty(dt))
    dt = .01;
  end
  
  if (~exist('alph')||isempty(alph))
    alph= .1;
  end  
   
  if (nargin <= 5)
    t0 = 0;
  end  
   
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
        proj = [1 2 3];
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
     
  if (isfield(S,'X0plot_opt'))
    Xopt = S.X0plot_opt;
  else
    Xopt = {'.', 'MarkerSize',6};
  end

  Xopt = {Xopt{:}, 'Color',[0 0 0]};
  %S.traj_plot_opt = {'g','LineWidth',2}
  S.Xfplot_opt = {'or','MarkerSize',6};
  S.X0plot_opt = {'+g','MarkerSize',6};
  hold on;
  nb_traj= numel(S.traj);
  
  switch(numel(proj))          
   case 2
    for i = 1:nb_traj            
      
      traj = S.traj{i};            
            
      if (numel(dt)==1)
        time = t0:dt:traj.time(end);
        if (time(end)<traj.time(end))
          time = [time traj.time(end)];
        end
      else 
        time = dt;
      end  
      
      Expat = interp1(traj.time', traj.Expa', time')';
      Xt =  interp1(traj.time', traj.X', time')';
      
      for k = 1:numel(time) 
        
        DX = [Expat(proj,k)'];
        X = [Xt(proj,k)'];
        plot(X(1),X(2),Xopt{:});
        
        if (k==1)
          rect(X-DX,2*DX,'b',.1);
        else
          rect(X-DX,2*DX,c,alph);    
        end
        
      end  
    end
    
 %   SplotTraj(S);

   otherwise
    for i = 1:nb_traj
      traj = S.traj{i};
      
      if (numel(dt)==1)
        time = t0:dt:traj.time(end);
        if (time(end)<traj.time(end))
          time = [time traj.time(end)];
        end
      else 
        time = dt;
      end  
    
      Expat = interp1(traj.time', traj.Expa', time')';
      Xt =  interp1(traj.time', traj.X', time')';
      
      
      for k = 1:numel(time)
        DX = Expat(proj,k)';
        X = Xt(proj,k)';
        plot3(X(1),X(2),X(3),Xopt{:});
        voxel(X-DX,2*DX,'b',.1);
      end  
      
    end
    
%    SplotTraj(S);

  end
  
