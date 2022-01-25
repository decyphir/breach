classdef sinusoid_signal_gen < signal_gen
    % sinusoid_signal_gen   A class derived from signal_gen to generate simple sinusoid signals.
    %  
    % spike_signal_gen Methods
    %   spike_signal_gen -  constructor, takes signal names, an optional signal interpolation 
    %                       methods for each signals or all signals, and an optional p0.
    %                       Each signal 'x' gets a 'x_sin_base', 'x_sin_amp',
    %                       'x_sin_decay', and 'x_sin_freq' parameter, with default 
    %                       values set as 0, 1, 1, and 0.2 respectively.
    %           
    %                         
    %  See also signal_gen.  
   
    methods 
        
        function this = sinusoid_signal_gen(signals)
           this.signals = signals; 
           this.params = {};
           for i_s = 1:numel(this.signals)
               this.params = { this.params{:} [this.signals{i_s} '_sin_base'] ...
                              [this.signals{i_s} '_sin_amp']... 
                              [this.signals{i_s} '_sin_decay']... 
                              [this.signals{i_s} '_sin_freq']};
               this.p0(4*(i_s-1)+1:4*i_s, 1) = [0 1 0 2];
           end

           this.params_domain = repmat(BreachDomain(), 1, numel(this.params));
           this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));     
        end
      
            
        function [X, time] = computeSignals(this, p, time) % compute the signals
            if numel(p) ~= 4*numel(this.signals)
                error('Wrong number of parameters for computing constant signal.' )
            end

            if size(p,1) ==1
                p = p';
            end

            for i_s = 0:numel(this.signals)-1 
                % with the variable order: base, amp, decay, freq                   
                pi_s = p(4*i_s+1:4*i_s+4);         
                X(i_s+1,:) = pi_s(1) + pi_s(2) * exp(pi_s(3)*time) .* sin(pi_s(4)*time);                
            end
           
        end
        
        function type = getType(this)
            type = 'sinusoid';
        end
    end
            
end
