classdef FalsificationProblem < BreachProblem
    %  FalsificationProblem A class dedicated to falsification of STL formulas.
    %
    %  FalsificationProblem Properties
    %    BrSet_False -  BreachSet updated with falsifying parameter vectors
    %                   and traces whenever some are found
    %    X_false     -  parameter values found falsifying the formula
    %    StopAtFalse - (default: true) if true, will stop as soon as a falsifying
    %                   parameter is found.
    %
    %  FalsificationProblem Methods
    %    GetBrSet_False - returns BrSet_False
    %
    % See also BreachProblem
    
    properties
        BrSet_False
        X_false
        obj_false
        StopAtFalse=true
    end
    
    methods
        
        % Constructor calls parent constructor
        function this = FalsificationProblem(BrSys, phi, params, ranges)
           
            Br = BrSys.copy();
            switch nargin
                case 2
                    params = Br.GetSysVariables();
                    req_params = Br.GetReqVariables();
                    if ~isempty(req_params)
                      Br.ResetDomain(req_params);
                    end
                    super_args{1} = Br;
                    super_args{2} = phi;
                    if isempty(params)
                        error('FalsificationProblem:invalid_system_variables', 'No valid system or input variables.');
                    end
                    super_args{3} = params;
                    
                case 3
                    super_args{1} = Br;
                    super_args{2} = phi;
                    super_args{3} = params;
                case 4
                    super_args{1} = Br;
                    super_args{2} = phi;
                    super_args{3} = params;
                    super_args{4} = ranges;
            end
            this = this@BreachProblem(super_args{:});
            this.obj_best=inf;
        end
        
        function ResetObjective(this)
            ResetObjective@BreachProblem(this);
            this.X_false = [];
            this.BrSet_False = [];
            this.obj_best = inf;
        end
        
        function obj = objective_fn(this,x)
            % For falsification, default objective_fn is simply robust satisfaction of the least
            robs = this.robust_fn(x);
            NaN_idx = isnan(robs); % if rob is undefined, make it inf to ignore it
            robs(NaN_idx) = inf;
            obj = min(robs);
        end     
        % Nothing fancy - calls parent solve then returns falsifying params
        % if found.
        function [Xfalse, res] = solve(this)
            res = solve@BreachProblem(this);
            Xfalse = this.X_false;
        end
        
        % Logging
        function LogX(this, x, fval)
            
            % Logging default stuff
            this.LogX@BreachProblem(x, fval);
            
            %  Logging falsifying parameters found      
            [~, i_false] = find(fval<0);
            if ~isempty(i_false)
                this.X_false = [this.X_false x(:,i_false)];                              
                if (this.log_traces)
                    if isempty(this.BrSet_False)
                        this.BrSet_False = this.BrSys.copy();
                    else
                        this.BrSet_False.Concat(this.BrSys);
                    end
                end
            end    
        end
        
        function b = stopping(this)
            b =  (this.time_spent >= this.max_time) ||...
                (this.nb_obj_eval>= this.max_obj_eval) ||...
                (this.StopAtFalse&&this.obj_best<0);
        end
        
        function BrFalse = GetBrSet_False(this)
            BrFalse = [];
            if this.log_traces
                BrFalse = this.BrSet_False;
            else
                [~, i_false] = find(this.obj_log<0);
                if ~isempty(i_false)
                    BrFalse = this.BrSys.copy();
                    BrFalse.SetParam(this.params, this.X_log(:, i_false));
                    if ~isempty(this.BrSys.log_folder)
                        BrFalse.Sim();
                    end
                end
                
            end
        end
        
        
        function DispResultMsg(this)
            this.DispResultMsg@BreachProblem();
                    
            if ~isempty(this.X_false)
                fprintf('Falsified with obj=%g\n', this.obj_best(end));
            else
                fprintf('No falsifying trace found.\n');
            end
        end
        
    end
end