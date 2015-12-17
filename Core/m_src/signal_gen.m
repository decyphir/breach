classdef signal_gen
    
    % Abstract class to generate a multidimensional signal of some type 
    properties
        signals % names of the signals it generates
        params  % parameters such as control points, etc
        p0      % default values  
    end
    
    
    methods
        function params = getParamNames(this)
            params = this.params;
        end
    end

    methods (Abstract)
        computeSignals(this,p, time) % compute the signals
    end

    
end

