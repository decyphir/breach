classdef req_monitor < signal_gen
    properties
        name
        R  % BreachRequirement object 
    end
    
    methods 
        function this =  req_monitor(name, signals, params, p0)
            if nargin ==0 
                this.name = 'req_null';                
            else
                this.name = name;
                if exist('params', 'var')
                    this.params = params;
                end
                
                if exist('signals', 'var')
                    this.signals_in = signals;
                end
                
                if exist('p0', 'var')
                    this.p0 = p0;
                end
            end            
        end
        
        function [t, Xout] = computeSignals(t, Xin, pin)
            Xout=[];
        end
        
        function status = set_mode(this, flag1, flag2)                 
            status = 0;
        end
        
        function [v,t,X] = eval(this, t,X,p)
            v =0; 
        end
    end    
    
end