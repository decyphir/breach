classdef expr_output_gen < stl_monitor
    
    methods
        function this = expr_output_gen(name, expr)
            formula = STL_Formula(name, [expr '>0']);
            this = this@stl_monitor(formula);
            this.signals ={name};
        end
    
         function [time, Xout] = computeSignals(this, time, X, p)
            this.init_tXp(time,X,p);                      
           [time, Xout] = this.get_standard_rob(this.formula, time);
         end
    
    
    end

    
end