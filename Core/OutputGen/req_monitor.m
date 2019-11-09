classdef req_monitor < output_gen
    properties
        name
    end
    
    methods 
        function [t, Xout] = computeSignals(t, Xin, pin)
            Xout=[];
        end
        
        function status = set_mode(this, flag1, flag2)                 
            status = 0;
        end
                
    end
    methods (Abstract)
         eval(this, t, X,p)
    end
    
end