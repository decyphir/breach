classdef pulse_signal_gen < signal_gen
    
    % Pulse signal generation 
    
    methods 
        % default pulse has base value 0, period 1, duty width 0.5 and amp
        % 1.
        function this = pulse_signal_gen(signals)
           this.signals = signals; 
           this.params = {};
           this.p0 = zeros( 4*numel(signals), 1 );
           for i_s = 1:numel(this.signals)
               this.params = { this.params{:} [this.signals{i_s} '_base_value'] ...
                              [this.signals{i_s} '_pulse_period']... 
                              [this.signals{i_s} '_pulse_width']...
                              [this.signals{i_s} '_pulse_amp']};
               this.p0(4*(i_s-1)+1:4*i_s,1) = [0 1 0.5 1];          
           end     
        end
                            
        function X = computeSignals(this,p, time) % compute the signals
            nb_p = 4; 
            if numel(p) ~= nb_p*numel(this.signals)
                error('Wrong number of parameters for computing constant signal.' )
            end
            if size(p,1) ==1
                p = p';
            end
            
            X = repmat(p(1:nb_p:end), 1, numel(time));
            
            unshifted_pulse = @(T,d,t) sin(((2*pi)/T)*t) > cos(pi*d);
            advance = @(T,d) T*(0.25 - 0.5*d);
            pulse = @(T,d,t) unshifted_pulse(T,d, t + advance(T,d));
            
            for i_s = 0:numel(this.signals)-1 
                pi_s = p(nb_p*i_s+1:nb_p*i_s+nb_p);
                base = pi_s(1);
                period = pi_s(2);
                duty = pi_s(3);
                amp = pi_s(4); 
                X(i_s+1,:) = base + amp*pulse(period, duty, time);              
            end
           
        end
        
        function type = getType(this)
            type = 'pulse';
        end
    end
            
end


