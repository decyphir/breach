classdef param_constraint_monitor < stl_monitor  % maybe should be the other way around ..? 
    
    methods
        function this= param_constraint_monitor(varargin)
        % param_constraint_monitor(phi) assumes that phi does not contain
        % temporal operators and signal
            this = this@stl_monitor(varargin{:});
            this.signals = {};  % no output signal
        end
        
        
        function [tau, Xout] = computeSignals(this, t, X, p, tau)
            % does nothing
            tau = t;
            Xout = X;
        end
        
        function [v, Xout] = eval(this, t, X,p)
           if nargout >1
               Xout = X;
           end
           this.assign_params(p);
           expr = get_expr(this.formula);
           v = eval(expr); 
        
        end
        
    end
end