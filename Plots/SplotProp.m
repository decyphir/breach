function [val, XI, YI, ZI] = SplotProp(Pf, prop, opt)
% SPLOTPROP  Plots the satisfaction of a prop wrt param values in a
% param set. Works for one (1d function) or two  (2d surface) varying  parameters. 
%  
%  Synopsys:  [val, XI, YI, ZI] = SplotProp(Pf, prop, opt)  
%
%  Prerequisite: prop was evaluated before
%   
%  Returns the plotted values.
%   
  
%% read options 
 
  if (~exist('opt','var'))
    opt = [];
  end
  
  if isfield(opt, 'style')
    style = opt.style;
    if isempty(style)
       style =  {'-b'};
    end
  else
    style = {'-b'};   
  end
  
  if isfield(opt, 'contour')  
    use_contour = opt.contour;
  else
    use_contour = 0;
  end
  
  if isfield(opt, 'nb_contour')  
    nb_contour = opt.nb_contour;
  else
    nb_contour = 0;
  end
  
  
%% plot the thing
  iprop = find_prop(get_id(prop), Pf.props_names);  
  
  if iprop
    switch numel(Pf.dim)
     case 1
      val = cat(1, Pf.props_values(iprop,:).val);
      val = val(:,1);
      [x, order] = sort(Pf.pts(Pf.dim,:));
      val = val(order);
      plot(x,val, style{:});
      xlabel(Pf.ParamList{Pf.dim},'Interpreter','none');  
      ylabel(disp(prop, -1),'Interpreter','none');
     
     case 2
      val = cat(1, Pf.props_values(iprop,:).val);
      if isfield(opt,'plot_pts')
        plot3(Pf.pts(Pf.dim(1),:), Pf.pts(Pf.dim(2),:), val+.01*abs(val), 'MarkerSize', 14);
      else
          
        
        Z = val(:,1);        
        [XI, YI, ZI] = QuickMeshSf(Pf,Z);
        
        if use_contour 
          clf;
          contourf(XI,YI,ZI, 256, 'EdgeColor', 'None'); 
          hold on;
          if (max(val)*min(val)<0)
            contour(XI,YI,ZI, [0 0], 'LineWidth', 2, 'EdgeColor', 'k');
            legend('sat/unsat', 'boundary') ;
          end
                    
        end
          
        xlabel(Pf.ParamList{Pf.dim(1)},'Interpreter','none');  
        ylabel(Pf.ParamList{Pf.dim(2)},'Interpreter','none');
        zlabel('Quantitative Satisfaction');
        title(disp(prop, -1),'Interpreter','none');                        

        
        % 0 level-set contour
        hold on;
        [c h] = contour(XI,YI,ZI,  [0 0], 'LineWidth',2,'LineColor','k');
        clabel(c,h);
        
        % Redraw for colormap
        prop_cmap(val);
        caxis([min(val) max(val)]);
        colorbar;
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
 
  
  
  
