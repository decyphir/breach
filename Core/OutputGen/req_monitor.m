classdef req_monitor < output_gen
    properties
        name
    end
    
    methods 
        function [t, Xout] = computeSignals(t, Xin, pin)
            Xout=[];
        end
    end
    methods (Abstract)
         eval(this, t, X,p)
    end
    
end