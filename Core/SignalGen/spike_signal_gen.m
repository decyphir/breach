classdef spike_signal_gen < signal_gen
    % spike_signal_gen   A class derived from signal_gen to generate simple spike signals.
    %  
    % spike_signal_gen Methods
    %   spike_signal_gen -  constructor, takes signal names, an optional signal interpolation 
    %                       methods for each signals or all signals, and an optional p0.
    %                       Each signal 'x' gets a 'x_spike_base', 'x_spike_amp',
    %                       'x_spike_width', and 'x_spike_time' parameter, with default 
    %                       values set as 0, 1, 1, and 1 respectively. For spike_amptitude,
    %                       positive value result in upward spike while negative value  
    %                       represents downward spike. 
    %           
    %                         
    %  See also signal_gen.  
 
    properties
        interp_method;
    end
    methods 
        
        function this = spike_signal_gen(signals, method, p0)
           if ischar(signals)
              signals= {signals}; 
           end
            
           this.signals = signals;
           this.params = {};
           this.interp_method = {};

           if exist('method','var')
              if ischar(method) 
                  method = {method};
              end
           else
               method = repmat( {'spline'}, 1, numel(signals));
           end
           
           if ~exist('p0', 'var')
               p0=[];
           end
           
           num_m = length(method);
           if num_m ~= 1 || length(method) ~= numel(this.signals)
               error(['Interpolation method can either be set for each individual' ...
                   'signal or for all signals']);
           end
           this.interp_method = method;
      
           for i_s = 1:numel(this.signals)
               this.params = { this.params{:} [this.signals{i_s} '_spike_base'] ...
                              [this.signals{i_s} '_spike_amp']... 
                              [this.signals{i_s} '_spike_width']... 
                              [this.signals{i_s} '_spike_time']};
               this.p0(4*(i_s-1)+1:4*i_s, 1) = [0 1 1 1];
           end
           
           if ~isempty(p0)
               if length(p0)==length(this.p0)
                   this.p0(:) = p0(:);
               else
                   error('spike_signal_gen:bad_p0','Incorrect dimensions for p0');
               end
           end
           

           this.params_domain = repmat(BreachDomain(), 1, numel(this.params));
           this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));     
        end
      
            
        function [X, time] = computeSignals(this, p, time) % compute the signals
            if numel(p) ~= 4*numel(this.signals)
                error('Wrong number of parameters for computing spike signal.' )
            end

            if size(p,1) ==1
                p = p';
            end
            
            % set all signals as a constant signal of spike_base value
            X = repmat(p(1:4:end), 1, numel(time));
            method = this.interp_method{1};
            for i_s = 0:numel(this.signals)-1 
                % with the variable order: base, amp, width, time
                if length(this.interp_method) ~= 1
                    method = this.interp_method{i_s+1};
                end
                    
                pi_s = p(4*i_s+1:4*i_s+4);
                
              %  if time(1) > pi_s(4) - pi_s(3)/2.0 || time(end) < pi_s(4) + pi_s(3)/2.0
              %      warning('Spike is outside the signal range.' )
              %  end 

                i_left = find(time>=(pi_s(4) - pi_s(3)/2.0), 1);
                if  time(end) < pi_s(4) + pi_s(3)/2.0
                    i_right = length(time);
                else 
                    i_right = find(time>=(pi_s(4) + pi_s(3)/2.0), 1);
                end
                
                x = [pi_s(4) - pi_s(3)/2.0; pi_s(4); pi_s(4) + pi_s(3)/2.0];
                y = [0; pi_s(2); 0];
                x_i = time(i_left: i_right)';
                y_i = interp1(x,y,x_i,method);              
                X(i_s+1, i_left:i_right) = X(i_s+1, i_left: i_right)+ y_i';

            end
           
        end
        
        function type = getType(this)
            type = 'spike';
        end
        function args = getSignalGenArgs(this)
            args = {'interp_method'};         
        end
        
    end
            
end


