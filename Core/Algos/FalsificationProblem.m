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
        Rio_Mode
        Rio_Mode_log=[]
        val_max=inf
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
         
            if nargin>=1
                Br = BrSys.copy();
            end
            switch nargin
                case 0
                    super_args = {};
                case 2
                    [params, ipr] = Br.GetSysVariables();
                    params = params(ipr>Br.P.DimX);
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
        
            
        function set_IO_robustness_mode(this, mode, cap)
            switch mode
                case 'default'
                    this.robust_fn = @(x) (this.Spec.Eval(this.BrSys, this.params, x));
                    this.constraints_fn = [];
                case 'random'
                    this.robust_fn = @(x) this.boolean_verdict(x);
                    this.constraints_fn = [];
                case 'in'
                    this.robust_fn = @(x) (this.Spec.Eval_IO('in', 'rel', this.BrSys, this.params, x));
                    if exist('cap','var')
                        this.val_max = cap;
                    end
                    this.constraints_fn = [];
                case 'out'
                    this.robust_fn = @(x) (this.Spec.Eval_IO('out', 'rel', this.BrSys, this.params, x));
                    if exist('cap','var')
                        this.val_max = cap;
                    end
                    this.constraints_fn = [];
                case 'constrained'
                    this.robust_fn = @(x) (this.Spec.Eval_IO('out', 'rel', this.BrSys, this.params, x));
                    if exist('cap','var')
                        this.val_max = cap;
                    end
                    this.constraints_fn = @(x) (this.Spec.Eval_IO('in', 'abs', this.BrSys, this.params, x));
                case 'combined'
                    this.robust_fn = @(x) (this.combined_IO_robustness(x));
                    this.constraints_fn = [];
            end
        end
        
        function ert = boolean_verdict(this, x)
            rob = this.Spec.Eval(this.BrSys, this.params, x);
            if rob >= 0
                ert = 0.5;
            else
                ert = -0.5;
            end
        end 

        function rio = combined_IO_robustness(this, x)
            ri = this.Spec.Eval_IO('in', 'abs', this.BrSys, this.params, x);
            ro = this.Spec.Eval_IO('out','rel', this.BrSys, this.params, x);
            rio = atan(ri) + atan(ro);
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
                varname = this.whoamI;
                warning('off','MATLAB:Figure:FigureSavedToMATFile');
                if ~evalin('base', ['exist(''' varname ''', ''var'')'])
                    assignin('base', varname,this);
                    evalin('base', ['save(''' FileSave ''',''' varname ''');']);
                    evalin('base', ['clear ' varname ';']);
                else
                    evalin('base', ['save(''' FileSave ''',''' varname ''');']);
                end
                warning('on','MATLAB:Figure:FigureSavedToMATFile');                
            end
        end
        
        
        % Logging
        function LogX(this, x, fval, cval, x_stoch)
        % LogX  log variable parameter value tested by optimizers
       
            % Logging default stuff
            this.LogX@BreachProblem(x, fval, cval,x_stoch);
            % this.Rio_Mode_log = [ this.Rio_Mode_log this.Rio_Mode ];  % not sure what this is useful for
            
            %  Logging falsifying parameters found
            [~, i_false] = find(min(fval)<0);
            if ~isempty(i_false)
                this.X_false = [this.X_false x(:,i_false)];
                this.obj_false = [this.obj_false fval(:,i_false)];
                if (this.log_traces)%&&~this.use_parallel&&~(this.BrSet.UseDiskCaching)  % FIXME - logging flags and methods need be revised
                    if isempty(this.BrSet_False)
                        if ~isempty(this.Spec.BrSet)
                            this.BrSet_False = this.Spec.BrSet.copy();
                        end
                    else
                        this.BrSet_False.Concat(this.BrSys, true);
                    end
                end
            end
            
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
            if ~strcmp(this.display, 'off')
                
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
end