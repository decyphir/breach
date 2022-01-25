classdef exponential_signal_gen < signal_gen
    % exponential_signal_gen   A class derived from signal_gen to generate simple  
    %                          exponential signals.
    %  
    % Methods
    %   exponential_signal_gen -  constructor, takes signal names and an optional p0.
    %                       Each signal 'x' gets a 'x_exp_base', 'x_exp_amp',
    %                       and 'x_exp_rate' parameter, with default values set as
    %                       1, 10, and 0.1 respectively. Note that when
    %                       'x_exp_amp' value is positive, the class will
    %                       generate exponentially increase signals while
    %                       negative value will result in exponentially
    %                       decreasing signals.
    %                       
    %           
    %                         
    %  See also signal_gen.  

    methods 
        
        function this = exponential_signal_gen(signals,p0)
           
           if ischar(signals)
               signals = {signals};
           end
           this.signals = signals; 
           this.params = {};
      
           if ~exist('p0', 'var')
               p0=[];
           end
           
           for i_s = 1:numel(this.signals)
               this.params = { this.params{:} [this.signals{i_s} '_exp_base'] ...
                   [this.signals{i_s} '_exp_amp']...
                   [this.signals{i_s} '_exp_rate']};
                   this.p0(3*(i_s-1)+1:3*i_s, 1) = [1 10 0.1];
           end
           if ~isempty(p0)
               if length(p0)==length(this.p0)
                   this.p0(:) = p0(:);
               else
                   error('exponential_signal_gen:bad_p0','Incorrect dimensions for p0');
               end
           end
           
           this.params_domain = repmat(BreachDomain(), 1, numel(this.params));
           this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));     
        end
      
            
        function [X, time] = computeSignals(this, p, time) % compute the signals
            if numel(p) ~= 3*numel(this.signals)
                error('Wrong number of parameters for computing exponential signal.' )
            end

            if size(p,1) ==1
                p = p';
            end
            
            % set all signals as a constant signal of exp_base value
            %X = repmat(p(1:3:end), 1, numel(time));

            for i_s = 0:numel(this.signals)-1 
                % with the variable order: base, amp, rate                   
                pi_s = p(3*i_s+1:3*i_s+3);         
                X(i_s+1,:) = pi_s(1)+pi_s(2)*(1-exp(-pi_s(3).*time));                
            end
           
        end
        
        function type = getType(this)
            type = 'exponential increase';
        end
    end
            
end
