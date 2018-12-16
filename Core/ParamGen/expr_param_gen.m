classdef expr_param_gen < param_gen
    % expr_param_gen applies a generic expression 
    
    properties
        expr
    end
    methods
        function this = expr_param_gen(expr, param_out, p0)
            this.expr = expr;
            this.params = extract_expr_params(expr);
            this.params_out = {param_out};
            
            if nargin<3
                p0 = 0;
            end
            this.p0 = p0;
            
        end
        
        function p_out = computeParams(this, p_in)
            p_out = NaN(numel(this.params_out), size(p_in, 2));
            for ip = 1:size(p_in,2)
                this.assign_params(p_in(:,ip));
                p_out(:,ip) = eval(this.expr);
            end
        end
        
    end
    
    
    
end