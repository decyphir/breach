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
        StopAtTrue=1
        obj_true
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
        % here obj is the opposite of robustness, so we still try to make it negative. 

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
            objs = -robs; % possible array of values corresponding to multiple objectives
            NaN_idx = isnan(objs); % if undefined, make it -inf to ignore it
            objs(NaN_idx) = -inf;
            obj = max(objs,[],1)'; % we want to make all of them negative 
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
            this.obj_true = [];
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

            % Logging default stuff
            this.LogX@BreachProblem(x, fval, cval,x_stoch);
            % this.Rio_Mode_log = [ this.Rio_Mode_log this.Rio_Mode ];  % not sure what this is useful for

            %  Logging satisfying parameters found
            [~, i_true] = find(max(fval)<0);  % for max sat, we need all fval to be negative, so all robustness are positive 
            if ~isempty(i_true)
                this.X_true = [this.X_true x(:,i_true)];
                this.obj_true = [this.obj_true fval(:,i_true)];
                if (this.log_traces)%&&~this.use_parallel&&~(this.BrSet.UseDiskCaching)  % FIXME - logging flags and methods need be revised
                    if isempty(this.BrSet_True)
                        if ~isempty(this.Spec.BrSet)
                            this.BrSet_True = this.Spec.BrSet.copy();
                        end
                    else
                        this.BrSet_True.Concat(this.BrSys, true);
                    end
                end
            end

        end

        function b = stopping(this)
           
            b =  this.stopping@BreachProblem();

            if this.StopAtTrue~=0 && ... % StopAtTrue is not false
                    rem(this.StopAtTrue,1)==0 &&...  % StopAtTrue is integer
                    this.StopAtTrue<=sum(this.obj_log<0) % found enough
                b = true;
            end

        end
        
        
        
        function [BrTrue, BrTrue_Err, BrTrue_badU] = GetBrSet_True(this)
            BrTrue = [];
            if this.log_traces&&~this.use_parallel 
                BrTrue = this.BrSet_True;
            else
                [~, i_true] = find(max(fval,[],1)<0);  % for max sat, we need all fval to be negative, so all robustness are positive
                if ~isempty(i_true)
                    BrTrue = this.BrSys.copy();
                    BrTrue.SetParam(this.params, this.X_log(:, i_true));
                    if this.BrSys.UseDiskCaching
                        BrTrue.Sim();
                    end
                end

                [BrTrue, BrTrue_Err, BrTrue_badU] = this.ExportBrSet(BrTrue);

            end
                        
        end
        
       
        function DispResultMsg(this)
            fprintf('\n ---- Best value %g found at\n', this.obj_best);
            param_values = this.x_best;
            for ip = 1:numel(this.params)
                fprintf( '        %s = %g\n', this.params{ip},param_values(ip))
            end
            fprintf('\n');
        end
        
        
    end
   % Should we display -fval or +fval ? For now, for consistency with falsif, we display fval, which we always minimize 
   % methods (Access=protected)
        %function display_status(this, fval, cval, obj_best)
        %  if ~strcmp(this.display,'off')
        %        if nargin==1
        %            obj = this.Spec.val;
        %            if ~isempty(this.Spec.precond_monitors)
        %                cval = min(min(this.Spec.traces_vals_precond)); % bof bof
        %            else
        %                cval = NaN;
        %            end
        %        end
        %       if nargin<4
        %           obj_best= -this.obj_best;
        %       end
        %
        %       
        %        display_status@BreachProblem(this, obj, cval, obj_best);
        %   end
        %end
    %end
end