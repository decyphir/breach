function [val XI YI ZI] = SplotProp(Pf, prop,opt)
% SPLOTPROP  Plots the quantitative satisfaction of property prop for the different
%  param values in Pf. Works in 1d and 2d only. 
%  
%  Synopsys:  val = SplotProp(Pf, prop, opt)  
%
%  Prerequisite: prop was evaluated before
%   
%  Returns the plotted values.
%   
  
  if (~exist('opt','var'))
    opt.style = {'-b'};
  end
    
  iprop = find_prop(get_id(prop), Pf.props_names);  
  
  if iprop
    switch numel(Pf.dim)
     case 1
      val = cat(1, Pf.props_values(iprop,:).val);
      val = val(:,1);
      plot(Pf.pts(Pf.dim,:),val,opt.style{:});
      xlabel(Pf.ParamList{Pf.dim},'Interpreter','none');  
      ylabel(disp(prop, -1),'Interpreter','none');
     
     
     case 2
      val = cat(1, Pf.props_values(iprop,:).val);
      if isfield(opt,'plot_pts')
        plot3(Pf.pts(Pf.dim(1),:), Pf.pts(Pf.dim(2),:), val+.01*abs(val), 'MarkerSize', 14);
      else
        Z = val(:,1);

        [XI YI ZI] = QuickMeshSf(Pf,Z);         
        
        xlabel(Pf.ParamList{Pf.dim(1)},'Interpreter','none');  
        ylabel(Pf.ParamList{Pf.dim(2)},'Interpreter','none');
        zlabel(disp(prop, -1),'Interpreter','none');        
      end      
      
     otherwise 
      return;
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
 
  
  