classdef step_signal_gen < signal_gen
    % step_signal_gen   A class derived from signal_gen to generate simple step signals.
    %  
    % step_signal_gen Methods
    %   step_signal_gen -  constructor, takes signal names, and an optional p0.
    %                       Each signal 'x' gets a 'x_base_value', 'x_step_time',
    %                       and 'x_step_amp' parameter, with default value
    %                       to 0, 1, and 1 respectively. 
    %                         
    %  See also signal_gen.  
 
    
    methods 
        
        function this = step_signal_gen(signals, p0)
           this.signals = signals; 
           this.params = {};
           for i_s = 1:numel(this.signals)
               this.params = { this.params{:} [this.signals{i_s} '_base_value'] ...
                              [this.signals{i_s} '_step_time']... 
                              [this.signals{i_s} '_step_amp']};
               this.p0(3*(i_s-1)+1:3*i_s,1) = [0 1 1];

           end     
           
           if nargin == 2
               this.p0 = p0;
           end

           this.params_domain = repmat(BreachDomain(), 1, numel(this.params));
           this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));
 
           
        end
        
                
        function params = getParamNames(this) % get parameterization names, e.g., signal1_u0, signal2_u0, etc                                             
            params = this.params;
        end
            
        function [X, time] = computeSignals(this,p, time) % compute the signals
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


