classdef expr_output_gen < stl_monitor
    
    methods
        function this = expr_output_gen(name, expr)
            formula = STL_Formula(name, [expr '>0']);
            this = this@stl_monitor(formula);
        end
    end
    
end