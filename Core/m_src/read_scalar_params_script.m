function [params___, p0___] = read_scalar_params_script(script___)
% [params, p0] = read_scalar_params_script() reads scalar variables in the
% workspace and their values


run(script___);
list_var___ = who;

params___ = {};
p0___ = [];

for iv___ = 1:numel(list_var___)
    
    var___= list_var___{iv___};
    
    b = eval(sprintf('isscalar(%s)&&isnumeric(%s)' ,var___,var___));
    if b
        params___ = [params___ {var___}];
        p0___     = [p0___ eval(var___)];
    end

end
    

end
