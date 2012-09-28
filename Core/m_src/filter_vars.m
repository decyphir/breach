function [vars vals] =  filter_vars(mdl, exclude)
%
% FILTER_VAR filter variables found in Simulink models : exclude capitalized
%   constants and change lookup tables into single variables
%
% Syntax: [vars vals] =  filter_var(VarsOut, exclude)
%
%   exclude is a regular expression patterns that should not be in the
%   names of the variables
%  
%   Ex: [vars vals] = filter_vars( 'model', '[A-Z]') will exclude all
%   variable with capitalized letters in them
%
%

  MAX_EL = 100;
  VarsOut = Simulink.findVars(mdl, 'WorkspaceType', 'base');
  newp = 0; % counts parameters found

  if nargin == 1
    exclude= 'tspan';
  end
  
  vars = {'dumb_par_ignore'};
  vals = [0];
  
  for i = 1:numel(VarsOut)
    var  = VarsOut(i);     
    vname = var.Name;
    if ~isempty(regexp(vname, exclude))
%      fprintf('excluding %s\n', vname);
      continue;
    end
    v = evalin('base',vname);
    if (isscalar(v))      
      newp= newp+1;
      vars{newp}= vname;
      vals(newp) = v;            
    else
      fprintf('found %s, %d %d table\n', vname,size(v,1), size(v,2) );         
      if numel(v) > MAX_EL
        if input(['Wait, this has ' num2str( numel(v)) ' elements, are you ' ...
                            ' sure to continue (0 or 1) ?'])
          for  i = 1:size(v,1)
            for j= 1:size(v,2)
              if (v(i,j)~=0)   
                newp=newp+1;
                vars{newp} = [vname '_tab_' num2str(i) '_' num2str(j)];
                vals(newp) = v(i,j);                    
              end
            end
          end   
          
        end
      end
    end                 
  end
