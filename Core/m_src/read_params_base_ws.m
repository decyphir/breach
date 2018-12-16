function [params, p0] = read_params_base_ws()
% [params, p0] = read_params_base_ws() reads scalar variables in the
% workspace and their values

exclude = {'ans'}; % what else? 
list_var = setdiff(evalin('base', 'who'), exclude);
params = {};
p0 = [];

for iv = 1:numel(list_var)
    
    var= list_var{iv};
    b = evalin('base',sprintf('isscalar(%s)&&isnumeric(%s)' ,var,var));
    if b
        params = [params {var}];
        p0 = [p0 evalin('base', var)];
    end

end
    

end