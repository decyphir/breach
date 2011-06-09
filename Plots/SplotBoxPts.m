function SplotBoxPts(S,proj, ipts, opt, col,alph)
%
%   Plots the points in the field pts of a parameter set and the boxes they represent.
%
%   SplotBoxPts(S, proj, ipts, opt, col,alpha)
%    
%   Inputs: 
%   
%    -  S        Box sampling set. Any set with a field 'Pts'
%                will do though.
%
%    -  proj     chooses the parameters to plot; can be numbers or
%                parameters names; all if [] 
%
%    -  ipts       indices of the pts to plot; all if absent or []   
% 
%    -  opt      Uses the plotting options defined in field X0plot_opt or
%                default if this is absent 
%    -  col      Color (e.g. 'g' 'r' 'k' etc)
%
%    -  alph     alpha rendering
%
%   Outputs:     
%      
%    -  Some figure, hopefully an interesting one.
% 
  
  
  if (~exist('proj')||isempty(proj))
    proj = 1:min(3, numel(S.dim));
    proj = S.dim(proj);
  end
  
  
  if (~exist('ipts')||isempty(ipts))
    ipts=1:size(S.pts,2);  
  end  
  
  if (exist('opt')&&(~isempty(opt)))
    SplotPts(S,proj, ipts, opt);
  else
    SplotPts(S,proj, ipts);
  end
  
  nb_pts = numel(S.pts(1,:));
  if (nb_pts==0)
    return;
  end
    
  if (~exist('col'))
    col = 'b';
  end
  
  if (~exist('alph'))
    alph = .1;
  end

   
  hold on;

  if (isfield(S,'plot_proj'))
    proj = S.plot_proj;
  else
    if (~exist('proj')||isempty(proj))
      switch numel((S.dim))
       case {1}
        proj=S.dim(1);
       case {2}
        proj=S.dim(1 :2);
       otherwise
        proj =S.dim(1:3);
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
  
  switch(numel(proj))
      
   case 1 
    
    hold on;    
    if isfield(S,'ParamList')            
      xlabel(S.ParamList{proj(1)});  
    else
      xlabel(['x_' num2str(proj(1))]);
    end
    
    DX(2) = 0;
    i = find(S.dim == proj);
    
    if (~isempty(i))          
      
      for k=ipts
      
        DX(1) = S.epsi(i,k);
        X = [S.pts(proj,k)'-DX(1),-DX(2)];               
        rect(X,2*DX,col,alph);
        
      end
    end
%    set(gca, 'YLim', [-1 1], 'YtickLabel', {});
    
   case 2
    hold on;
    nb_pts= size(S.pts,2);
   
    for k=ipts      
      for j = 1:2
        i = find(S.dim == proj(j));
        if isempty(i)          
          DX(j) = 0;
        else
          DX(j) = S.epsi(i,k);
        end
      end        
      X = S.pts(proj,k)'-DX;
      rect(X,2*DX,col,alph);
    end
    
   case 3    
    hold on
    if isfield(S,'ParamList')            
      xlabel(S.ParamList{proj(1)});  
      ylabel(S.ParamList{proj(2)});
      zlabel(S.ParamList{proj(3)});
    else
      xlabel(['x_' num2str(proj(1))]);
      ylabel(['x_' num2str(proj(2))]);
      zlabel(['x_' num2str(proj(3))]);
    end
    nb_pts= size(S.pts,2);

    for k=ipts
      for j = 1:3
        i = find(S.dim == proj(j));
        if (isempty(i)||S.epsi(i,k)==0)          
          DX(j) = 0;         
        else
          DX(j) = S.epsi(i,k);
        end
      end
 
      X = S.pts(proj,k)'-DX;
      voxel(X,2*DX,col,alph);
    end
  end
  grid on;
  hold off;