classdef step_signal_gen
    
    % Class to single step multidimensional signals 
    properties 
        signals  % names of the signals it generates
        params % parameters such as number of control points, etc
    end
    
    methods 
        
        function this = step_signal_gen(signals)
           this.signals = signals; 
           this.params = {};
           for i_s = 1:numel(this.signals)
               this.params = { this.params{:} [this.signals{i_s} '_base_value'] ...
                              [this.signals{i_s} '_step_time']... 
                              [this.signals{i_s} '_step_amp']};
            end     
        end
                
        function params = getParamNames(this) % get parameterization names, e.g., signal1_u0, signal2_u0, etc                                             
            params = this.params;
        end
            
        function X = computeSignals(this,p, time) % compute the signals
            if numel(p) ~= 3*numel(this.signals)
                error('Wrong number of parameters for computing constant signal.' )
            end
            if size(p,1) ==1
                p = p';
            end
            
            X = repmat(p(1:3:end), 1, numel(time));
            
            for i_s = 0:numel(this.signals)-1 
                pi_s = p(3*i_s+1:3*i_s+3);
                i_after = find(time>pi_s(2));
                X(i_s+1,i_after) = X(i_s+1,i_after)+pi_s(3);              
            end
           
        end
        
        function type = getType(this)
            type = 'step';
        end
    end
            
end


