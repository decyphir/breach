classdef pulse_signal_gen < signal_gen
    % pulse_signal_gen   A class derived from signal_gen to generate periodic pulse signals.
    %
    % pulse_signal_gen Methods
    %   pulse_signal_gen -  constructor, takes signal names, and an optional p0.
    %                       Each signal 'x' gets a 'x_base_value', 'x_pulse_period',
    %                       'x_pulse_width', 'x_pulse_amp' and 'x_pulse_delay' parameter, with default
    %                       to 0, 1, 0.5, 1 and 0, respectively.
    %
    %  See also signal_gen.
    
    methods
        function this = pulse_signal_gen(signals, p0)
            if ischar(signals)
                signals = {signals};
            end
            this.signals = signals;
            this.params = {};
            this.p0 = zeros( 5*numel(signals), 1 );
            for i_s = 1:numel(this.signals)
                this.params = { this.params{:} [this.signals{i_s} '_base_value'] ...
                    [this.signals{i_s} '_pulse_period']...
                    [this.signals{i_s} '_pulse_width']...
                    [this.signals{i_s} '_pulse_amp']...
                    [this.signals{i_s} '_pulse_delay']};
                this.p0(5*(i_s-1)+1:5*i_s,1) = [0 1 0.5 1 0];
            end
            
            if nargin==2
                this.p0 =p0;
            end
            
           this.params_domain = repmat(BreachDomain(), 1, numel(this.params));
           this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));
 
            
        end
        
        function [X, time] = computeSignals(this,p, time) % compute the signals
            nb_p = 5;
            if numel(p) ~= nb_p*numel(this.signals)
                error('Wrong number of parameters for computing constant signal.' )
            end
            if size(p,1) ==1
                p = p';
            end
            
            X = repmat(p(1:nb_p:end), 1, numel(time));
            
            unshifted_pulse = @(T,d,t) sin(((2*pi)/T)*(t+eps)) > cos(pi*d);
            advance = @(T,d) T*(0.25 - 0.5*d);
            pulse = @(T,d,t) unshifted_pulse(T,d, t + advance(T,d));
            
            for i_s = 0:numel(this.signals)-1
                pi_s = p(nb_p*i_s+1:nb_p*i_s+nb_p);
                base = pi_s(1);
                period = pi_s(2);
                duty = pi_s(3);
                amp = pi_s(4);
                delay = pi_s(5);
                if delay ==0
                    X(i_s+1,:) = base + amp*pulse(period, duty, time);
                elseif delay>0
                    time_not_0 = time(time>= delay);
                    X(i_s+1, time>=delay) = base + amp*pulse(period, duty, time_not_0-time_not_0(1));
                else 
                    dt = time(2)-time(1);
                    time_shift = [time(1)+delay:dt:time(1)-dt time];
                    pulse_s =  base + amp*pulse(period, duty, time_shift-delay);
                    X(i_s+1,:) = pulse_s(time_shift>=time(1));                
                end
                
            end
            
            function type = getType(this)
                type = 'pulse';
            end
        end
        
    end
end  
    
