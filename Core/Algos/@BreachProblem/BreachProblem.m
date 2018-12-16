classdef BreachProblem < BreachStatus
    % BreachProblem A class for generic optimization problem involving STL specifications
    %
    % A BreachProblem is essentially created from a BreachSystem and a
    % property. E.g.,
    %
    %   problem = BreachProblem(BrSys, phi);
    %
    % In that case BrSys parameter ranges determines the variables and
    % search domains. BrSys parameter vectors are used as initial values
    % for the problem. Alternatively, parameters to optimize and ranges can
    % be specified explictly using the syntax:
    %
    %   problem = BreachProblem(BrSys, phi, params, ranges);
    %
    % where params is a cell array of parameters and ranges a corresponding
    % array of ranges.
    %
    % When a problem is created, the objective function is constructed from
    % the robust satisfaction of the property phi. A solver can then be
    % selected to solve it via the solve method.
    %
    % BreachProblem Properties (inputs)
    %   BrSet          -  BreachSet used as (initial) domain for optimization
    %   BrSys          -  BreachSystem used by the solver to compute new
    %                     traces or new satifaction values
    %   solver         -  default: 'basic', use list_solvers to get a list of available solvers
    %   solver_options -  option structure for the solver. See each solver
    %                     help to know available options (E.g., for matlab solvers, this is
    %                     often set using the optimset command).
    %   max_time       -  maximum wall-time budget allocated to optimization
    %   max_obj_eval   -  maximum number of objective function evaluation (
    %   often translates into number of simulations of system under test)
    %   log_traces     -  (default=true) logs all traces computed during optimization                 
    %   T_spec         -  time 
    % 
    % BreachProblem Properties (outputs)
    %   BrSet_Logged    -  BreachSet with parameter vectors used during optimization.
    %   BrSet_Best      -  BreachSet with the best parameter vector found during optimization.
    %   res             -  a result structure, specific to each solver
    %   X_log, obj_log  -  all values tried by the solver and corresponding objective function values.
    %   x_best,obj_best -  best value found
    %
    %
    % BreachProblem Methods
    %   list_solvers    - returns a list of strings of available solvers
    %   setup_solver    - setup the solver given as argument with default options
    %                     and returns these options.
    %   solve           - calls the solver 
    %   GetLog          - returns a BreachRequirement object with logged
    %   traces
    %   GetBest         - returns a BreachRequirement object with the best
    %   trace (worst satisfaction in case of falsification) 
    %
    % See also FalsificationProblem, ParamSynthProblem, ReqMiningProblem
    
    %% Properties 
    properties
        objective
        x0
        solver= 'global_nelder_mead'   % default solver name
        solver_options    % solver options
        Spec              % BreachRequirement object reset to R0 for each objective evaluation
        T_Spec=0
        
        constraints_fn    % constraints function
        robust_fn         % base robustness function - typically the robust satisfaction of some property by some trace
    end
    
    % properties related to the function to minimize
    properties
        BrSet
        BrSys         % BreachSystem reset for each objective evaluation
        BrSet_Best   
        BrSet_Logged
        R0            % BreachRequirement object initial 
        R_log    % BreachRequirement object logging requirement evaluations during solving 
        
        params
        domains
        lb
        ub
        Aineq
        bineq
        Aeq
        beq
        res
        X_log
        obj_log
        x_best
        obj_best   = inf
        log_traces = true
    end
    
    % misc options
    properties
        display = 'on'
        freq_update = 10 
        use_parallel = 0
        max_time = inf
        time_start = tic
        time_spent = 0
        nb_obj_eval = 0
        max_obj_eval = 100
        mixed_integer_optim_solvers = {'ga'};
    end
    
    %% Static Methods
    methods (Static)
       
        function solvers = list_solvers(verbose)
        % list_solvers display the list of (supposedly) available solvers.      
        % TODO check dependency on locally installed toolboxes 
        
        solvers = ...
            {'init', ...
            'basic',...
            'global_nelder_mead',...
            'binsearch',...
            'fminsearch',...
            'cmaes'...
            };
        
        solvers_others = ...
            {'fmincon', ...
            'simulannealbnd', ...
            'optimtool',...
            'ga',...
            };
        
        if ~exist('verbose')||(verbose~=0)
            for i_solv = 1:numel(solvers)
                if isequal(solvers{i_solv}, 'global_nelder_mead')
                disp([ solvers{i_solv} ' (default)']);
                else
                    disp(solvers{i_solv});
                end
            end
        end
        
        if ~exist('verbose')||(verbose~=0)
            for i_solv = 1:numel(solvers_others)
                if exist(solvers_others{i_solv})
                    disp(solvers_others{i_solv});
                    solvers= [solvers solvers_others{i_solv}];
                end
            end
        end
        
        end
    end
    
    %% Methods
    methods
             
        %% Constructor
        function this = BreachProblem(BrSet, phi, params, ranges)
            
            if ~isa(phi, 'BreachRequirement')
                phi = BreachRequirement(phi);
            end
            
            this.R0 = phi.copy(); 
            this.Spec = phi;
            
            switch nargin 
                case 2
                    this.ResetObjective(BrSet);
                case 3
                    this.ResetObjective(BrSet, params);
                case 4
                    this.ResetObjective(BrSet, params, ranges);
            end
            % setup default solver
            this.setup_solver();
            
            % reset display
            rfprintf_reset();
            
        end
        
        function Reset_x0(this)
            
            x0__ = zeros(numel(this.params), size(this.BrSet.P.pts,2));
            for ip = 1:numel(this.params)
                x0__ip =  this.BrSet.GetParam(this.params{ip});
                if isempty(x0__ip)
                    x0__ip =  this.Spec.GetParam(this.params{ip});
                end
                if isempty(x0__ip)
                    error('BreachProblem:unknown_param', ['Parameter ' this.params{ip} ' is neither a system parameter nor a requirement parameter.']);
                end
                x0__(ip,:) = x0__ip;
            end
            
            this.BrSet.SetParam(this.params, x0__,'spec');  % not sure this is useful anymore, if ever
            this.x0 = unique(x0__', 'rows')';% remove duplicates, I guess. 
            
        end
        
        function ResetObjective(this, BrSet, params, ranges)
            if nargin == 1
                BrSet = this.BrSet;
            else
                this.BrSet = BrSet.copy();
            end
            
            this.BrSet.Sys.Verbose=0;
                  
            % Parameter ranges
            if ~exist('params','var')
                params = this.BrSet.GetBoundedDomains();
            else
                if ischar(params)
                    params = {params};
                end
            end
            this.params= params;
            this.domains = BrSet.GetDomain(params);
            
            if ~exist('ranges', 'var')
                ranges = BrSet.GetParamRanges(params);
                lb__ = ranges(:,1);
                ub__ = ranges(:,2);
            else
                lb__ = ranges(:,1);
                ub__ = ranges(:,2);
                this.BrSet.SetParam(params, 0.5*(ranges(:,2)-ranges(:,1)), true); % adds parameters if they don't exist 
                this.BrSet.ResetDomains();
                this.BrSet.SetDomain(params, 'double', ranges);
            end
           
            this.lb = lb__;
            this.ub = ub__;            
            
            this.Reset_x0;
            
            % Reset display
            rfprintf_reset();
            
            % robustness
            this.BrSys = BrSet.copy();
            
            this.BrSys.SetParam(this.params, this.x0(:,1), 'spec');
            [~, ia] = unique( this.BrSys.P.pts','rows');
            if numel(ia)<size(this.BrSys.P.pts,2)
                this.BrSys = this.BrSys.ExtractSubset(ia);
            end

            this.BrSys.Sys.Verbose=0;
            
            this.robust_fn = @(x) (this.Spec.Eval(this.BrSys, this.params, x));
            
            % objective function
            this.objective = @(x) (objective_wrapper(this,x));      
            
            this.BrSet_Best = [];
            this.BrSet_Logged = [];
            this.res = [];
            this.X_log = [];
            this.obj_log= [];
            this.obj_best =inf;
            this.time_spent = 0;
            this.nb_obj_eval = 0;
        end
           
        function ResetTimeSpent(this)
            this.time_spent= 0; 
            this.time_start = tic; 
        end
           
        %% Options for various solvers
        function [solver_opt, is_gui] = setup_solver(this, solver_name, is_gui)
            if ~exist('solver_name','var')
                solver_name = this.solver;
            end
            if exist('is_gui', 'var')&&is_gui==true
                try
                    solver_opt = eval(['this.setup_' solver_name '(true);' ]);
                catch
                    is_gui = false;
                    solver_opt = eval(['this.setup_' solver_name '();' ]);
                end
            else
                solver_opt = eval(['this.setup_' solver_name '();' ]);
                is_gui = false;
            end
        end
        
        function solver_opt = setup_init(this)
            solver_opt = struct();
            this.solver= 'init';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_optimtool(this)
            solver_opt = optimset('Display', 'iter');
            this.display = 'off';
            solver_opt.lb = this.lb;
            solver_opt.ub = this.ub;
            this.solver = 'optimtool';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_fmincon(this)
            disp('Setting options for fmincon solver');
            solver_opt = optimoptions('fmincon', 'Display', 'iter');
            this.display = 'off';
            if this.max_obj_eval < inf
                solver_opt = optimoptions(solver_opt, 'MaxFunEvals', this.max_obj_eval);
            end
            if this.use_parallel 
                solver_opt = optimoptions(solver_opt,'UseParallel', true);
                if ~this.BrSys.use_parallel
                    this.BrSys.SetupParallel();
                end
            end
            this.solver = 'fmincon';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_fminsearch(this)
            disp('Setting options for fminsearch solver');
            solver_opt = optimset('fminsearch');
            solver_opt = optimset(solver_opt, 'Display', 'iter');
            if this.max_obj_eval < inf
                solver_opt = optimset(solver_opt, 'MaxFunEvals', this.max_obj_eval);
            end
            if this.max_time < inf
                solver_opt = optimset(solver_opt, 'MaxTime', this.max_time);
            end
            solver_opt = optimset(solver_opt, 'MaxIter', numel(this.params)*200);
            this.display = 'off';
            % Mathworks currently don't support parallel fminsearch 
            this.solver = 'fminsearch';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_simulannealbnd(this)
            disp('Setting options for simulannealbnd solver');
            this.display = 'off';
            solver_opt = saoptimset('Display', 'iter');
            if this.max_time < inf
                solver_opt = saoptimset(solver_opt, 'MaxTime', this.max_time);
            end
            if this.max_obj_eval < inf
                solver_opt = saoptimset(solver_opt, 'MaxFunEvals', this.max_obj_eval);
            end
            % Mathworks currently don't support parallel simulannealbnd
            
            this.solver = 'simulannealbnd';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_cmaes(this)
            disp('Setting options for cmaes solver - use help cmaes for details');
            solver_opt = cmaes();
            solver_opt.Seed = 0;
            solver_opt.LBounds = this.lb;
            solver_opt.UBounds = this.ub;
            if this.max_obj_eval < inf
                solver_opt.MaxFunEvals = this.max_obj_eval;
            end
            this.display = 'off';
            solver_opt.SaveVariables = 'off';
            solver_opt.LogModulo = 0;
            %solver_opt.DispModulo = 0; % need to disable when running
            %multiple cmaes instances
            if this.use_parallel 
                solver_opt.EvalParallel = 'yes';
            end
            this.solver = 'cmaes';
            this.solver_options = solver_opt;
        end
        
        %% solve functions for various solvers
        function res = solve(this)
            
            % reset display
            rfprintf_reset();
            
            % reset time
            this.ResetTimeSpent();
            
            % create problem structure
            problem = this.get_problem();
            
            switch this.solver
                case 'init'
                    this.display_status_header();
                    res = FevalInit(this);
                    
                case 'basic'
                    res = this.solve_basic();
                    
                case 'global_nelder_mead'
                    res = this.solve_global_nelder_mead();
                    
                case 'cmaes'
                    % adds a few more initial conditions
                    nb_more = 10*numel(this.params)- size(this.x0, 2);
                    if nb_more>inf % what is this for? Not sure
                        Px0 = CreateParamSet(this.BrSet.P, this.params,  [this.lb this.ub]);
                        Px0 = QuasiRefine(Px0, nb_more);
                        this.x0 = [this.x0' GetParam(Px0,this.params)]';
                    end
                    
                    [x, fval, counteval, stopflag, out, bestever] = cmaes(this.objective, this.x0', [], this.solver_options);
                    res = struct('x',x, 'fval',fval, 'counteval', counteval,  'stopflag', stopflag, 'out', out, 'bestever', bestever);
                    
                case 'ga'
                    res = solve_ga(this, problem);
                    
                case 'fmincon'
                    while ~this.stopping
                        [x,fval,exitflag,output] = feval(this.solver, problem);
                        res = struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);
                        if ~this.stopping % restart
                            problem.x0 = this.generate_new_x0;
                        end
                    end
                    
                case 'fminsearch'
                    while ~this.stopping
                        if this.use_parallel
                            num_works = this.BrSys.Sys.Parallel;
                            for idx = 1:num_works
                                problem.x0 = this.generate_new_x0;
                                F(idx) = parfeval(this.solver, 4, problem);
                            end
                            res = cell(1,num_works);
                            for idx = 1:num_works
                                [completedIdx, x, fval,exitflag,output] = fetchNext(F);
                                this.nb_obj_eval = this.nb_obj_eval + output.funcCount;
                                res{completedIdx} = struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);
                            end
                        else
                            [x,fval,exitflag,output] = feval(this.solver, problem);
                            res = struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);
                            if ~this.stopping % restart
                                problem.x0 = this.generate_new_x0;
                            end
                        end
                    end
                case 'simulannealbnd'
                    if this.use_parallel
                        num_works = this.BrSys.Sys.Parallel;
                        for idx = 1:num_works
                            F(idx) = parfeval(this.solver, 4, problem);
                        end
                        res = cell(1,num_works);
                        for idx = 1:num_works
                            [completedIdx, x, fval,exitflag,output] = fetchNext(F);
                            this.nb_obj_eval = this.nb_obj_eval + output.funccount;
                            res{completedIdx} = struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);
                        end
                    else
                        [x,fval,exitflag,output] = feval(this.solver, problem);
                        res = struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);
                    end
                    
                case 'optimtool'
                    problem.solver = 'fmincon';
                    optimtool(problem);
                    res = [];
                    return;

                case 'binsearch'
                    res = solve_binsearch(this);

                otherwise
                    res = feval(this.solver, problem);
                    
            end
            this.res = res;
            this.DispResultMsg(); 
        
            %% Saving run in cache folder
            this.SaveInCache();
        
        end
        
        function SaveInCache(this)
            if this.BrSys.UseDiskCaching
                FileSave = [this.BrSys.DiskCachingRoot filesep class(this) '_Runs.mat'];
                evalin('base', ['save(''' FileSave ''',''' this.whoamI ''');']);
            end
        end
        
        %% Utility functions for solvers
        
        % function res = FevalInit(this,X0)
        % defined in external file
        
        function X0 = init_basic_X0(this)
            % returns initial vectors
            BrQ = this.BrSet.copy();
            BrQ.ResetParamSet();
            BrQ.SetParamRanges(this.params, [this.lb this.ub])
            BrC = BrQ.copy();
            nb_corners = this.solver_options.nb_max_corners;
            nb_samples = this.solver_options.nb_new_trials;
            step = this.solver_options.start_at_trial;
            
            BrC.P = CreateParamSet(BrC.Sys,this.params,[this.lb this.ub]);
            BrC.CornerSample(nb_corners);
            XC = BrC.GetParam(this.params);
            nb_corners= size(XC, 2);
            qstep = step-nb_corners;
            if qstep>=0
                % skips corners
                BrQ.QuasiRandomSample(nb_samples, step);
                X0 = BrQ.GetParam(this.params);
            else
                qnb_samples = nb_samples+qstep;
                if qnb_samples>0  % needs to finish corners plus some
                    BrQ.QuasiRandomSample(qnb_samples);
                    XQ = BrQ.GetParam(this.params);
                    X0 = [XC(:,step+1:end) XQ];
                else % more corners than samples anyway
                    X0 = XC(:,step+1:end);
                end
                
            end
            
        end
        
        function x0 = generate_new_x0(this)
            x0 = (this.ub-this.lb).*rand(length(this.ub),1) + this.lb;
        end
        
        function problem = get_problem(this)
            problem = struct('objective', this.objective, ...
                'fitnessfcn', this.objective, ... % for ga
                'x0', this.x0, ...   
                'nvars', size(this.x0, 1),... % for ga
                'solver', this.solver,...
                'Aineq', this.Aineq,...
                'bineq', this.bineq,...
                'Aeq', this.Aeq,...
                'beq', this.beq,...
                'lb', this.lb,...
                'ub', this.ub,...
                'nonlcon', this.constraints_fn,...
                'intcon',[],...
                'rngstate',[],...
                'options', this.solver_options);
            
            % Checks whether some variables are integer
            for ip = 1:numel(this.params)
                dom = this.BrSys.GetDomain(this.params{ip});
                if strcmp(dom.type, 'int')  
                    problem.intcon = [problem.intcon ip];
                end
            end
            
        end
        
        % check the MIP options for suppoted solvers and change the this.params
        function setup_mixed_int_optim(this, method)
            if ~exist('method')||isempty(method)
                method = 'map_enum_to_int';
            end
           
            switch method
                case 'map_enum_to_int'
                    
                    enum_idx = find(this.params_type_idx('enum'));
                    for ii = 1:length(enum_idx)
                        enum_param = this.params{enum_idx(ii)};
                        enum_domain = this.BrSet.GetDomain(enum_param);
                        enum_br_domain = BreachDomain('enum', enum_domain.enum);
                        pg = enum_idx_param_gen(enum_param, enum_br_domain);
                        this.BrSet.SetParamGen({pg});
                        this.params{enum_idx(ii)} = pg.params{1};
                        fprintf('Mapped %s to %s\n', pg.params_out{1}, pg.params{1} );
                    end
                    this.ResetObjective();
                case 'exhaustive'
                    idx = union(find(this.params_type_idx('enum')),find(this.params_type_idx('int')));
                    BrSet = this.BrSet.copy();
                    BrSet.SampleDomain(this.params(idx), 'all');
                    this.ResetObjective(BrSet, this.params(setdiff(1:numel(this.params), idx)));
                    fprintf('Sampled enum and int domains exhaustively, resulting in %g values for each objective evaluation.\n', BrSet.GetNbParamVectors());
            end
            
        end
                
        function idx = params_type_idx(this, ptype)
            domains = this.BrSys.GetDomain(this.params);
            idx = cell2mat(arrayfun(@(x) strcmp(x.type, ptype), ...
                domains, 'UniformOutput', false));
        end
       
        function add_constraint(this, phi)
            this.constraints_fn = @(x) (deal(-this.BrSys.GetRobustSat(phi, this.params, x, this.T_Spec), []));
        end
        
        %% Parallel 
        function SetupParallel(this, varargin)
            this.use_parallel = 1;  
            % Create parallel pool and get number of workers
            this.BrSys.SetupParallel(varargin{:});      
            
            % Enable DiskCaching
            this.SetupDiskCaching();

            % Possible need to change the optimization option
            % this.setup_solver();
            % TODO: review when needed on solver-by-solver basis - maybe
            % warning in order ?
            
            
        end
        
        function StopParallel(this)
            this.use_parallel = 0;
            this.BrSys.StopParallel();
        end
        
        function SetupDiskCaching(this, varargin)
            this.log_traces = 0;  % FIXME
            this.BrSys.SetupDiskCaching(varargin{:});
        end
        
        %% Objective wrapper        
        function obj = objective_fn(this,x)
            % default objective_fn is simply robust satisfaction 
            obj = this.robust_fn(x);
        end
        
        function fval = objective_wrapper(this,x)
             % objective_wrapper calls the objective function and wraps some bookkeeping           
             
             if size(x,1) ~= numel(this.params)
                x = x';
             end
        
            nb_eval =  size(x,2);
            fval = inf*ones(size(this.Spec.req_monitors,2), nb_eval);
            fun = @(isample) this.objective_fn(x(:, isample));
            nb_iter = min(nb_eval, this.max_obj_eval);
     
            if this.stopping()==false
                if nb_iter == 1 || ~this.use_parallel
                    for iter = 1:nb_iter

                        % checks whether x has already been computed or not
                        % skips logging if already computed 
                        % might cause inconsistencies down the road...
                        
                        if ~isempty(this.X_log)
                            xi = x(:, iter);
                            idx = find(sum(abs(this.X_log-repmat(xi, 1, size(this.X_log, 2)))) == 0,1);
                        else
                            idx=[];
                        end
                        if ~isempty(idx)
                            fval(:,iter) = this.obj_log(:,idx);
                        else
                            % calling actual objective function
                            fval(:,iter) = fun(iter);
                        
                            % logging and updating best
                            this.LogX(x(:, iter), fval(:,iter));
                            
                            % update status
                            if rem(this.nb_obj_eval,this.freq_update)==0
                                this.display_status();
                            end
                        end
                        % stops if falsified or other
                        if this.stopping()
                            break
                        end
                        
                    end
                else % Parallel case 
                    
                    % Launch tasks
                    for iter = 1:nb_iter
                        par_f(:,iter) = parfeval(fun,1, iter);
                    end
                    
                    fq = this.freq_update;                  
                    for iter=1:nb_iter
                        [idx, value] = fetchNext(par_f);
                        fval(:,idx) = value;
                        this.LogX(x(:, idx), fval(:,idx));
                        
                        % update status
                        if rem(iter,fq)==0
                            this.display_status();
                        end
                        if this.stopping()
                            cancel(par_f);
                            break
                        end
                    end
                end
            else
                fval = this.obj_best*ones(1, nb_eval);
            end
            
        end
        
        function b = stopping(this)
            b =  (this.time_spent >= this.max_time) ||...
                    (this.nb_obj_eval>= this.max_obj_eval);
        end
              
        %% Misc methods
        function x = CheckinDomain(this,x)
          for ip = 1:numel(this.params)
                x(ip) = this.domains(ip).checkin(x(ip));  
          end
        end
        
        function LogX(this, x, fval)
            % LogX logs values tried by the optimizer

            this.X_log = [this.X_log x];
            this.obj_log = [this.obj_log fval];
            
            if (this.log_traces)&&~(this.use_parallel)&&~(this.BrSet.UseDiskCaching) % FIXME - logging flags and methods need be revised
                if isempty(this.BrSet_Logged)
                    this.BrSet_Logged = this.BrSys.copy();
                    this.R_log = this.Spec.copy();
                else
                    this.BrSet_Logged.Concat(this.BrSys);
                    this.R_log.Concat(this.Spec);
                end
                this.Spec = this.R0.copy();
            end
            
            [fmin , imin] = min(min(fval));
            x_min =x(:, imin);
            if fmin < this.obj_best
                this.x_best = x_min;
                this.obj_best = fval(:,imin);
                this.BrSet_Best = this.BrSys.copy(); % could be more efficient - also suspicious when several x... 
            end
            
            % Timing and num_eval       
            this.nb_obj_eval= numel(this.obj_log(1,:));
            this.time_spent = toc(this.time_start);
           
        end
        
        function DispResultMsg(this)
       
            % display status of last eval if not already done
            if rem(this.nb_obj_eval,this.freq_update)
                this.display_status();
            end
            
            % DispResultMsg message displayed at the end of optimization
            if this.time_spent > this.max_time
                fprintf('\n Stopped after max_time was reached.\n');
            end
            
            if this.nb_obj_eval > this.max_obj_eval
                fprintf('\n Stopped after max_obj_eval was reached (maximum number of objective function evaluation.\n' );
            end
            
            if numel(this.res) > 1 && this.use_parallel
                fprintf('\nReports from different parallel optimization runs.\n');
                for idx = 1:numel(this.res)
                    fprintf('Run %d\n', idx);
                    this.Display_Best_Results(this.res{idx}.fval, this.res{idx}.x);
                    if this.obj_best > this.res{idx}.fval
                        this.obj_best = this.res{idx}.fval;
                        this.x_best = this.res{idx}.x; 
                    end
                end
                fprintf('\nIn summary, the best results among all runs: \n');
            end

            this.Display_Best_Results(this.obj_best, this.x_best);
            
        end
        
        function Display_Best_Results(this, best_fval, param_values)
            fprintf('\n ---- Best value %g found with\n', min(best_fval));
            
            for ip = 1:numel(this.params)
                % check the enum_idx variable and restore the params
                value = param_values(ip);
                fprintf( '        %s = %g\n', this.params{ip},value)
            end
            fprintf('\n');
        end
             
        function [BrOut, Berr, BbadU] = GetBrSet_Logged(this)
            % GetBrSet_Logged gets BreachSet object containing parameters and traces computed during optimization
            BrOut = this.BrSet_Logged;
            if isempty(BrOut)
                BrOut = this.BrSys.copy();
                BrOut.ResetSimulations();
                BrOut.SetParam(this.params, this.X_log);
                BrOut.Sim();
            end
            [BrOut, Berr, BbadU] = this.ExportBrSet(BrOut); 
            
        end
                
        function BrBest = GetBrSet_Best(this)
            BrBest = this.BrSet_Best;
            if isempty(BrBest)
                BrBest = this.BrSys.copy();
                BrBest.SetParam(this.params, this.x_best, 'spec');
                BrBest.Sim();
            end
            
            BrBest =  this.ExportBrSet(BrBest);
            
        end
              
         function summary = SaveResults(this, varargin)
            BLog = this.GetLog();   
            summary = BLog.SaveResults(varargin{:});
         end
        
        function Rlog = GetLog(this,varargin)
            Rlog = this.GetBrSet_Logged(varargin{:});
        end
        
        function Rbest = GetBest(this,varargin)
            Rbest = this.GetBrSet_Best(varargin{:});
        end
        
    end
    
    methods (Access=protected)
        
        function display_status_header(this)
            if ~isempty(this.Spec.precond_monitors)
                hd_st = sprintf(  '#calls (max:%5d)        time spent (max: %g)     [current  obj]     (current best)   [constraint]\n',...
                    this.max_obj_eval, this.max_time);
            else
                hd_st = sprintf(  '#calls (max:%5d)        time spent (max: %g)     [current  obj]     (current best) \n',...
                    this.max_obj_eval, this.max_time);
            end
            fprintf(hd_st);
        end
        
        function display_status(this,fval, const_val)
            
            if ~strcmp(this.display,'off')
                if nargin==1
                    fval = this.obj_log(:,end); % bof bof
                    if ~isempty(this.Spec.precond_monitors)
                        const_val = min(min(this.Spec.traces_vals_precond));
                    end
                end
                
                st__= sprintf('     %5d                    %7.1f               [%s]    (%s)', ...
                    this.nb_obj_eval, this.time_spent,  num2str(fval','%+5.5e '), num2str(this.obj_best', '%+5.5e '));
                if exist('const_val', 'var')
                    st__ = sprintf([st__ '       [%s]\n'], num2str(const_val', '%+5.5e '));
                else
                    st__ = [st__ '\n'];  
                end
                
                
                switch this.display
                    case 'on'
                        fprintf(st__);
                    case 'light'
                        rfprintf(st__);
                end
            end
        end
        
        function [BrOut, Berr,  BbadU] = ExportBrSet(this,B)
            % ExportBrSet prepares a BreachSet such as best, logged, etc to be
            % returned
            Berr = [];
            BbadU = [];
            if isempty(B)
                BrOut = [];
                return;
            end
             B.Sys.Verbose=1;
           
            [idx_ok, idx_sim_error, idx_invalid_input, st_status]  = B.GetTraceStatus();
            
            if ~isempty(idx_sim_error)||~isempty(idx_invalid_input)
                [B, Berr, BbadU] = FilterTraceStatus(B);
                this.disp_msg(['Warning: ' st_status],1);
            end
            
            if ~isempty(idx_ok)&&B.hasTraj();
                BrOut = this.Spec.copy();
                BrOut.ResetSimulations();
                BrOut.Eval(B);
            end
        end
                
        
    end
    
end
