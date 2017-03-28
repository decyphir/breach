classdef MaxSatProblem < BreachProblem
    % MaxSatProblem A variant of BreachProblem - maximize robustness and logs positive value
    %
    %  FalsificationProblem Properties
    %    BrSet_True -  BreachSet updated with falsifying parameter vectors
    %                   and traces whenever some are found
    %    X_True     -  parameter values found falsifying the formula
    %    StopAtTrue - (default: true) if true, will stop as soon as a falsifying
    %                   parameter is found.
    %
    %  FalsificationProblem Methods
    %    GetBrSet_True - returns BrSet_True
    %
    % See also BreachProblem
    
    properties
        BrSet_True
        X_True
        StopAtTrue=false
    end
    
    
    methods
        function this = MaxSatProblem(BrSys, phi, params, ranges)
            switch nargin
                case 2
                    super_args{1} = BrSys;
                    super_args{2} = phi;
                case 3
                    super_args{1} = BrSys;
                    super_args{2} = phi;
                    super_args{3} = params;
                case 4
                    super_args{1} = BrSys;
                    super_args{2} = phi;
                    super_args{3} = params;
                    super_args{4} = ranges;
            end
            this = this@BreachProblem(super_args{:});
        end
        
        function obj = objective_fn(this,x)
            obj = -min(this.robust_fn(x)); % maximizes the min robustness
        end
        
        function ResetObjective(this)
            ResetObjective@BreachProblem(this);
            this.X_True = [];
            this.BrSet_True = [];
        end
        
        
        % Nothing fancy - calls parent solve then display falsifying params
        % if found.
        function [XTrue, res] = solve(this)
            res = solve@BreachProblem(this);
            XTrue = this.X_True;
        end
        
        % Logging
        function LogX(this, x, fval)
            
            % Logging default stuff
            this.LogX@BreachProblem(x, fval);
            
            %  Logging satisfying parameters and traces
            if fval < 0
                this.X_True = [this.X_True x];
                if isempty(this.BrSet_True)
                    this.BrSet_True = this.BrSys.copy();
                else
                    this.BrSet_True.Concat(this.BrSys);
                end
                
                if this.StopAtTrue == true
                    this.stopping = true;
                end
            end
            
        end
        
        function BrTrue = GetBrSet_True(this)
            if this.log_traces
                BrTrue = this.BrSet_False;
            else
                [~, i_true] = find(this.obj_log>=0);
                if ~isempty(i_true)
                    BrTrue = this.BrSys.copy();
                    BrTrue.SetParam(this.params, this.X_log(:, i_true));
                    if ~isempty(this.BrSys.log_folder)
                        BrTrue.Sim();
                    end
                end
                
            end
        end
        
       
        function DispResultMsg(this)
            fprintf('\n ---- Best robustness value %g found at\n', -this.obj_best);
            param_values = this.x_best;
            for ip = 1:numel(this.params)
                fprintf( '        %s = %g\n', this.params{ip},param_values(ip))
            end
            fprintf('\n');
        end
        
    end
end