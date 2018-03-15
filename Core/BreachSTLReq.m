classdef BreachSTLReq < BreachConstraint
    
    methods
        function this = BreachSTLReq(formula)
            
            if nargin==0
                return
            else
                if isa(formula, 'char')||isa(formula, 'STL_Formula')
                    this.formula = stl_monitor(formula);
                elseif isa(formula, 'stl_monitor')
                    this.formula = formula;
                end
            end
        end
        
        function val = evalTrace(this,traj)
                [~,  val] = this.getRobustSignal(traj, 0); %  for each trace, we collect robustness at time 0
        end
        
        function [tau, val] = getRobustSignal(this,traj, tau)
                [tau,  val] = this.formula.computeSignals(traj.time,traj.X, traj.params, tau);
        end
        
        function disp(this)
            disp(disp(this.formula.formula));
        end
        
    end
    
end