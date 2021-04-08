classdef MyFalsificationProblem < FalsificationProblem
    % Example of a class deriving from FalsificationProblem specialized for one solver 
    methods
        function this = MyFalsificationProblem(BrSys, phi)
            this = this@FalsificationProblem(BrSys, phi); % calls parent constructor
            % add whatever is needed by your solver 
        end
        
        % setup options - will be called by the constructor
        function setup_solver(this)
            this.display = 'off';
            this.solver_options = optimset('Display', 'iter');
        end
        
        % custom objective function - robust_fn is obtained by the constructor from BrSys and phi
        function [fval, cval, x_stoch] = objective_fn(this,x)
            rob = min(this.robust_fn(x)); % note: robust_fn might return an array of values 
            fval = rob*norm(x);  % variation: we can try to maximize the norm of x            
            cval = inf;
            x_stoch = [];  % no stochastic variables
        end
        
        % call solver and return a structure with results. 
        function res = solve(this)
            this.ResetTimeSpent();   % set timer to 0
            [x, fval] = fmincon(this.objective, ... % objective calls objective_fn plus some bookkeeping
            this.x0, ...          % obtained from BrSys
             [], [], [],[], ...   % linear (in)-equalities
            this.lb, this.ub, ... % obtained by constructor from ranges 
            [], this.solver_options);
            res = struct('x_best',x, 'fval', fval); 
        end     
    end
end