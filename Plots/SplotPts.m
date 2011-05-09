function SplotPts(S,proj,ipts,opt)
%  
%  Plots parameters points in a parameter set S
%  
%  Usage:
%          SplotPts(S [, proj, ipts, opt])
%
%   Plots the points in field Pts      
%   
%   Inputs: 
%   
%    -  S        Box sampling set. 
%    -  proj     chooses the parameters to plot; can be numbers or
%                parameters names; all if [] 
%    -  ipts     indices of the pts to plot; all if absent or []   
%    -  opt      plot options; uses the plotting options defined in field
%                X0plot_opt or default if this is absent
%
%   see demoVDP1 for examples
  
 
  if (~exist('ipts')||isempty(ipts))
    ipts=1:size(S(1).pts,2);  
  end
  
  % dealing with plot options
  
  if (~exist('opt'))
    opt = [];
    if (isfield(S(1),'X0plot_opt'))
      opt.plot_opt = S(1).X0plot_opt;
    else
      opt.plot_opt = {'+m','MarkerSize',6};       
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
    opt.plot_opt=  {'+m','MarkerSize',6};       
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
      proj(i) = FindParam(S(1),stproj{i});
    end
  end  
  
  proj = proj(proj~=0);

  %  Deal with several Sampling sets
  
  if numel(S)>1
    hold on;
    nb = numel(S);
    colors = hsv(nb);
    for i=1:nb
      opt.plot_opt ={'+','MarkerSize',6,'Color',colors(i,:)};
      SplotPts(S(i),proj,[],opt);
    end    
    hold off
    return
  end
  
  switch (numel(proj))
    
   case {1}
    hold on;    
    if isfield(S,'ParamList')            
      xlabel(S.ParamList{proj(1)});  
    else
      xlabel(['x_' num2str(proj(1))]);
    end
      
    x = S.pts(proj(1),ipts)*rescale;
    plot(x,0*x,opt.plot_opt{:});
        
   case {2}
    hold on;  
   
    if isfield(S,'ParamList')            
      xlabel(S.ParamList{proj(1)});  
      ylabel(S.ParamList{proj(2)});
      
    else
      xlabel(['x_' num2str(proj(1))]);
      ylabel(['x_' num2str(proj(2))]);
      
    end
      
    x = S.pts(proj(1),ipts)*rescale;
    y = S.pts(proj(2),ipts)*rescale;
    plot(x,y,opt.plot_opt{:});  
    xlabel(S.ParamList{proj(1)});  
    ylabel(S.ParamList{proj(2)});   
   
  otherwise
    hold on;
    
    if isfield(S,'ParamList')            
      xlabel(S.ParamList{proj(1)});  
      ylabel(S.ParamList{proj(2)});
      zlabel(S.ParamList{proj(3)});
    else
      xlabel(['x_' num2str(proj(1))]);
      ylabel(['x_' num2str(proj(2))]);
      zlabel(['x_' num2str(proj(3))]);
    end
   
    x = S.pts(proj(1),ipts)*rescale;
    y = S.pts(proj(2),ipts)*rescale;
    z = S.pts(proj(3),ipts)*rescale;    
		plot3(x,y,z,opt.plot_opt{:});
      
  end
  grid on;
  hold off;
%  set(gca,'FontSize',14,'FontName','times');
 