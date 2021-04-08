classdef MaxSatProblem < BreachProblem
    % MaxSatProblem A variant of BreachProblem - maximize robustness and logs positive value
    %
    %  FalsificationProblem Properties
    %    BrSet_True -  BreachSet updated with falsifying parameter vectors
    %                   and traces whenever some are found
    %    X_true     -  parameter values found falsifying the formula
    %    StopAtTrue - (default: true) if true, will stop as soon as a falsifying
    %                   parameter is found.
    %
    %  FalsificationProblem Methods
    %    GetBrSet_True - returns BrSet_True
    %
    % See also BreachProblem
    
    properties
        BrSet_True
        X_true
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
        
        function [obj, cval, x_stoch] = objective_fn(this,x)

            % For falsification, default objective_fn is simply robust satisfaction of the least
            this.Spec.ResetEval();
            
            % checks stochastic domain
            if ~isempty(this.stochastic_params)
                this.BrStoch.ResetParamSet();
                this.BrStoch.SampleDomain(size(x, 2));
                x_stoch = this.BrStoch.GetParam(this.stochastic_params);
                this.BrSys.SetParam(this.stochastic_params, x_stoch);
            else
                x_stoch = nan(0, size(x,2)); 
            end
                        
            this.robust_fn(x);
            robs = this.Spec.traces_vals;
            
            if (~isempty(this.Spec.traces_vals_precond))
                precond_robs = robs;
                for itr = 1:size(this.Spec.traces_vals_precond,1)
                    precond_robs(itr) = min(this.Spec.traces_vals_precond(itr,:));
                    if  precond_robs(itr)<0
                        robs(itr,:)= -precond_robs(itr);
                    end
                end
            end
            robs = -robs;
            NaN_idx = isnan(robs); % if rob is undefined, make it inf to ignore it
            robs(NaN_idx) = -inf;
            obj = max(robs,[],1)';
            cval = inf;
            if (~isempty(this.Spec.traces_vals_precond))
                NaN_idx = isnan(precond_robs); % if rob is undefined, make it inf to ignore it
                precond_robs(NaN_idx) = inf;
                cval = min(precond_robs,[],1)';
            end

        end
        
        function ResetObjective(this, varargin)
            ResetObjective@BreachProblem(this, varargin{:});
            this.X_true = [];
            this.BrSet_True = [];
        end
        
        
        % Nothing fancy - calls parent solve then display falsifying params
        % if found.
        function [XTrue, res] = solve(this)
            res = solve@BreachProblem(this);
            XTrue = this.X_true;
        end
        
        % Logging
        function LogX(this, x, fval, cval, x_stoch)
            
            %  Logging satisfying parameters and traces
            [~, i_true] = find(fval>0);
            if ~isempty(i_true)
                this.X_true = [this.X_true x(:,i_true)];                              
                if (this.log_traces)&&~this.use_parallel
                    if isempty(this.BrSet_True)
                        this.BrSet_True = this.BrSys.copy();
                    else
                        this.BrSet_True.Concat(this.BrSys, true);
                    end
                end
            end
            
            % Logging default stuff
            this.LogX@BreachProblem(x, fval,cval, x_stoch);
            
        end
        
        function b = stopping(this)
            b =  this.stopping@BreachProblem();
            b= b||(this.StopAtTrue&&this.obj_best>0);
        end
        
        function [BrTrue, BrTrue_Err, BrTrue_badU] = GetBrSet_True(this)
            BrTrue = [];
            if this.log_traces&&~this.use_parallel 
                BrTrue = this.BrSet_True;
            else
                [~, i_false] = find(this.obj_log<0);
                if ~isempty(i_false)
                    BrTrue = this.BrSys.copy();
                    BrTrue.SetParam(this.params, this.X_log(:, i_false));
                    if this.BrSys.UseDiskCaching
                        BrTrue.Sim();
                    end
                end
                
                [BrTrue, BrTrue_Err, BrTrue_badU] = this.ExportBrSet(BrTrue);
                
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