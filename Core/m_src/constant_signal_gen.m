classdef constant_signal_gen < signal_gen
    
    % Class to generate multidimensional constant signals 
    
    methods         
        function this = constant_signal_gen(signals)
           this.signals = signals; 
           this.params = {};
           this.p0 = 0;
            for i_s = 1:numel(this.signals)
                this.params = { this.params{:} [this.signals{i_s} '_u0']};
            end     
        end
            
        function X = computeSignals(this,p, time) % compute the signals
            if numel(p) ~= numel(this.signals)
                error('Wrong number of parameters for computing constant signal.' )
            end
            if size(p,1) ==1
                p = p';
            end
            X = repmat(p, 1, numel(time)); 
        end
        
        function type = getType(this)
            type = 'constant';
        end
    end
            
end


