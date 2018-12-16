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
    
    methods (Static)
        function falsif_pb = load_runs(logFilePath)
            st = load([logFilePath, filesep, 'FalsificationProblem_Runs']);
            fn= fieldnames(st);
            falsif_pb = st.(fn{1});
            falsif_pb.SetupDiskCaching('DiskCachingRoot', logFilePath);
        end
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
        
        function ResetObjective(this, varargin)
            ResetObjective@BreachProblem(this,varargin{:});
            this.X_false = [];
            this.BrSet_False = [];
            this.obj_best = inf;
        end
        
        function obj = objective_fn(this,x)
            % For falsification, default objective_fn is simply robust satisfaction of the least
            this.robust_fn(x);
            robs = this.Spec.traces_vals;
            if (~isempty(this.Spec.traces_vals_precond))
                for itr = 1:size(this.Spec.traces_vals_precond,1)
                    precond_rob = min(this.Spec.traces_vals_precond(itr,:));
                    if  precond_rob<0
                        robs(itr,:)= -precond_rob;
                    end
                end
            end
            
            NaN_idx = isnan(robs); % if rob is undefined, make it inf to ignore it
            robs(NaN_idx) = inf;
            obj = min(robs,[],1)';
            
        end     
        
        % Nothing fancy - calls parent solve then returns falsifying params
        % if found.
        function [Xfalse, res] = solve(this)
            res = solve@BreachProblem(this);
            Xfalse = this.X_false;
        end
        
        function SaveInCache(this)
            if this.BrSys.UseDiskCaching
                FileSave = [this.BrSys.DiskCachingRoot filesep 'FalsificationProblem_Runs.mat'];
                evalin('base', ['save(''' FileSave ''',''' this.whoamI ''');']);
            end
        end
        
        
        % Logging
        function LogX(this, x, fval)
            %   LogX  log variable parameter value tested by optimizers
            
            %  Logging falsifying parameters found
            [~, i_false] = find(min(fval)<0);
            if ~isempty(i_false)
                this.X_false = [this.X_false x(:,i_false)];
                this.obj_false = [this.obj_false fval(:,i_false)];
                if (this.log_traces)&&~this.use_parallel&&~(this.BrSet.UseDiskCaching)  % FIXME - logging flags and methods need be revised
                    if isempty(this.BrSet_False)
                        this.BrSet_False = this.Spec.BrSet.copy();
                    else
                        this.BrSet_False.Concat(this.BrSys);
                    end
                end
            end
            % Logging default stuff
            this.LogX@BreachProblem(x, fval);
            
        end
        
        function b = stopping(this)
            b =  this.stopping@BreachProblem();
            b= b||(this.StopAtFalse&&any(this.obj_best<0));        
        end
        
        function [BrFalse, BrFalse_Err, BrFalse_badU] = GetFalse(this)
           BrFalse = this.BrSet_False;
            if isempty(BrFalse)
                [~, i_false] = find(min(this.obj_log)<0);
                if ~isempty(i_false)
                    BrFalse = this.BrSys.copy();
                    BrFalse.SetParam(this.params, this.X_log(:, i_false));
                    BrFalse.Sim();
                end
            end
            
            [BrFalse, BrFalse_Err, BrFalse_badU] = this.ExportBrSet(BrFalse);
        
        end
        
        function [BrFalse, BrFalse_Err, BrFalse_badU] = GetBrSet_False(this)
        % Use GetFalse. Keeping this for backward compatibility.
            [BrFalse, BrFalse_Err, BrFalse_badU] = this.GetFalse();
        end
        
        function DispResultMsg(this)
            this.DispResultMsg@BreachProblem();
            %if this.use_parallel && min(this.obj_best) < 0
            %    this.X_false = this.x_best;
            %end
            if ~isempty(this.X_false)
                fprintf('Falsified with obj = %g\n', min(this.obj_best(:,end)));
            else
                fprintf('No falsifying trace found.\n');
            end
        end
        
    end
end