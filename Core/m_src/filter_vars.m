function [vars, vals, ParamSrc] =  filter_vars(mdl, exclude)
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

load_system(mdl);
MAX_EL = 100;

%VarsOut = Simulink.findVars(mdl, 'WorkspaceType', 'base');
VarsOut = Simulink.findVars(mdl);
mdl_ws = get_param(mdl, 'modelworkspace');

newp = 0; % counts parameters found

if nargin == 1
    exclude= {'tspan'};
end

vars ={};
vals = [];
for i = 1:numel(VarsOut)
    var  = VarsOut(i);
    vname = var.Name;
    if any(~cellfun(@isempty, regexp(vname, exclude)))
        %      fprintf('excluding %s\n', vname);
        continue;
    end
    
    switch var.SourceType
        case 'base workspace'
            v = evalin('base',vname);
            v_src = 'base';
        case 'model workspace'
            v = getVariable(mdl_ws,var.Name);
            v_src = mdl_ws;
    end
    
    if (isnumeric(v))
        new_par(v, vname);
    elseif isstruct(v)
        fds = fieldnames(v);
        for svi = 1:numel(fds)
          new_par(v.(fds{svi}), [vname '.' fds{svi}]);
        end
    end
end

    function new_par(nv, nvname)
        if (isscalar(nv))
            newp= newp+1;
            vars{newp}= nvname;
            vals(newp) = nv;
        else
%         fprintf('found %s, %d %d table\n', vname,size(v,1), size(v,2) );  % TODO verbose option 
            if numel(nv) < MAX_EL
                for  i = 1:size(nv,1)
                    for j= 1:size(nv,2)
                        if (nv(i,j)~=0)
                            newp=newp+1;
                            vars{newp} = [nvname '_tab_' num2str(i) '_' num2str(j)];
                            vals(newp) = nv(i,j);
                            ParamSrc{newp} = v_src; % 'base' or model_ws 
                        end
                    end
                end
%            else
%              warning('Table has more than 100 element - ignored.') % TODO(?) force option
           end
        end
    end

end
