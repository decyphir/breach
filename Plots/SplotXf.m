function SplotXf(S,proj,ipts,opt)
%   Plots final points of computed trajectories in parameter set S
%
%   Usage:
%           SplotXf(S [, proj, ipts, opt])
%
%   Plots the points in field Xf (final points of computed trajectories)
%   in the current figure.
%      
%   
%   Inputs: 
%   
%    -  S        Parameter set. Any set with a field 'Xf'
%                will do though.
%    -  proj     chooses the parameters to plot; can be numbers or
%                parameters names; all if [] 
%    -  ipts       indices of the pts to plot; all if absent or []   
%    -  opt      Uses the plotting options defined in field Xfplot_opt or
%                default if this is absent 
%  

  
  if (~exist('ipts')||isempty(ipts))
    ipts=1:size(S(1).pts,2);  
  end
  

  if (~exist('opt'))
    opt = [];
    if (isfield(S(1),'X0plot_opt'))
      opt.plot_opt = S(1).X0plot_opt;
    else
      opt.plot_opt = {'+k','MarkerSize',6};       
    end
  else
    if (ischar(opt))    
      opt = {opt};
    end
  end
  
  if (iscell(opt))
    opt0 = opt;
    opt = [];      
    opt.plot_opt = opt0;
  end   

  % default plot option
  if (~isfield(opt, 'plot_opt'))
    opt.plot_opt=  {'+k','MarkerSize',6};       
  end    
  
  % rescaling axis
  
  rescale = 1;
  
  if isfield(opt, 'rescale')
    rescale = opt.rescale;
  end    
    
  if (isfield(S(1),'plot_proj'))
    proj = S(1).plot_proj;
  else
    if (~exist('proj')||isempty(proj))
      switch (numel(S(1).dim))
       case {1}
        proj=S(1).dim(1);
       case {2}
        proj=S(1).dim(1: 2);
       otherwise
        proj = S(1).dim(1:3);
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
  
     
  
  switch (numel(proj))
   case {1}
    hold on;
    x = S.Xf(proj(1),ipts);
    plot(x,0*x,opt.plot_opt{:});
    
   case {2}
    hold on;  
    x = S.Xf(proj(1),ipts);
    y = S.Xf(proj(2),ipts);
    plot(x,y,opt.plot_opt{:});
   
   otherwise
    hold on;
    x = S.Xf(proj(1),ipts);
    y = S.Xf(proj(2),ipts);
    z = S.Xf(proj(3),ipts);
    plot3(x,y,z,opt.plot_opt{:});            
    
  end

  grid on;