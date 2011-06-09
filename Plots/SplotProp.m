function val = SplotProp(Pf, prop,opt)

%  Plots the quantitative satisfaction of property prop for the different
%  param values in Pf. Works in 1d and 2d only. 
%  
%  Syntax:  val = SplotProp(Pf, prop, opt)  
%
%  Prerequisite: prop was evaluated before
%   
%  Returns the plotted values.
%   
  
  if (~exist('opt'))
    opt.style = {'-b'};
  end
  
  if (~isfield(opt,'rescale'))
    rescale = 1;
  else
    rescale = opt.rescale;
  end
  
  iprop = find_prop(get_id(prop), Pf.props_names);  
  
  if iprop
    switch numel(Pf.dim)
     case 1
      val = cat(1, Pf.props_values(iprop,:).val);
      val = val(:,1);      
      plot(Pf.pts(Pf.dim,:)*rescale,val,opt.style{:});
     case 2
      val = cat(1, Pf.props_values(iprop,:).val);
      if isfield(opt,'plot_pts')
        plot3(Pf.pts(Pf.dim(1),:), Pf.pts(Pf.dim(2),:), val+.01*abs(val), 'MarkerSize', 14);
      else
        Z = val(:,1);
        QuickMeshSf(Pf,Z);         
      
        if (rescale ~= 1)
          h = findobj(gca)
          for k =1:numel(h)
            try 
              X= get(h(k),'XData')*rescale;
              Y= get(h(k),'YData')*rescale;
              Z= get(h(k),'Zdata');
              break;
            end
          end
          cla;
          surf(X,Y,Z,'EdgeColor','None');        
        end
      end
    end
  end
  
  
function i = find_prop(st, props_names)
  i=0;
  for k = 1:numel(props_names)
    if strcmp(st,props_names{k})
      i = k;
      return;
    end    
  end
 
  
  