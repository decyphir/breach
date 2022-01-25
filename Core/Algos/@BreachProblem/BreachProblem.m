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
    %   max_obj_eval   -  maximum number of objective function evaluation 
    %                     (often translates into number of simulations of 
    %                      system under test)
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
    %                     traces
    %   GetBest         - returns a BreachRequirement object with the best
    %                     trace (worst satisfaction in case of falsification) 
    %
    % See also FalsificationProblem, ParamSynthProblem, ReqMiningProblem
    %
    
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
        
        %%Thao added
        insigma
    end
    
    % properties related to the function to minimize
    properties
        BrSet
        BrSys         % BreachSystem - reset for each objective evaluation
        BrSet_Best    
        BrSet_Logged   
        R0            % BreachRequirement object initial 
        R_log         % BreachRequirement object logging requirement evaluations during solving 
        
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
        
        % We can define the average robustness of the specification(s)
        % This is intended to be used for normaliazation
        % e.g. if we have two specs and the avg robustness values are 
        % [1 1000], we will in objective_wrapper divide the second specs
        % robustness by 1000 before passing it to the solver, in an effort
        % to keep the robustness values "equal"
        avgRobForNormalization
        
        % We define a "batch size" for parallel computations to not run out
        % of memory when doing parallel runs. If batch size is 50, then
        % solve() will calculate 50 objective function values in parallel,
        % then sort out the results, then calculate another 50 function
        % values in parallel etc. 
        parallelBatchSize = inf
        
        BrStoch       % Stochastic set of parameters, assigned to BrSys before every evaluation 
        stochastic_params      % list of stochastic parameters       
        stochastic_domains     % and corresponding domain
        X_stochastic_log       % logging stochastic values
        
        
    end
    
    % callbacks
    
    properties
       callback_obj  % called after each objective computation
       freq_update=1   % affects display and callbacks     
       is_paused = false
    end
    
    
    % misc options
    properties
        display = 'on'        
        use_parallel = 0
        max_time = inf
        time_start = tic
        time_spent = 0
        nb_obj_eval = 0
        max_obj_eval = 300
        num_constraints_failed = 0
        num_consecutive_constraints_failed = 0
        max_consecutive_constraints_failed = 100
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
            
            if nargin>0
                if ~isa(phi, 'BreachRequirement')
                    phi = BreachRequirement(phi);
                end
                
                % Initialize avg robustness to be used for normalization
                this.avgRobForNormalization = ones(numel(phi.req_monitors), 1);
                
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
            this.BrSet.verbose = 0; 
            
            % Parameter ranges
            if ~exist('params','var')
                [params, ipr] = this.BrSet.GetBoundedDomains();
                params = params(ipr>this.BrSet.P.DimX);
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
                this.BrSet.SetParam(params, 0.5*(ranges(:,2)+ranges(:,1)), true); % adds parameters if they don't exist 
                this.BrSet.ResetDomains(); % FIXME does not handle properly enum/int domains
                this.BrSet.SetDomain(params, 'double', ranges);
            end
           
            this.lb = lb__;
            this.ub = ub__;            
            
            this.Reset_x0;
            
            % Reset display
            rfprintf_reset();
            
            % robustness
            this.BrSys = this.BrSet.copy();
            
            this.BrSys.SetParam(this.params, this.x0(:,1), 'spec');
            [~, ia] = unique( this.BrSys.P.pts','rows');
            if numel(ia)<size(this.BrSys.P.pts,2)
                this.BrSys = this.BrSys.ExtractSubset(ia);
            end

            this.BrSys.Sys.Verbose=0;
            this.BrSys.verbose=0;
            
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
        
        
        %% Stochastic stuff
        function set_stochastic_params(this, sparams, sdomains)
            this.stochastic_params = sparams;
            this.stochastic_domains = sdomains;
            
            %[this.params, idx_non_stoch_params]= setdiff(this.params, sparams, 'stable'); 
            %this.domains = this.domains(idx_non_stoch_params);
            
            this.BrStoch = BreachSet(sparams);
            this.BrStoch.SetDomain(sparams,sdomains);
            
            % Reset objective function to consider only non stochastic 
            % parameters and their domains
            this.BrSet.SetEmptyDomain(sparams); 
            this.ResetObjective();            
            
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
        
        function solver_opt = setup_init(this,varargin)
            solver_opt = struct();
            this.solver= 'init';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_quasi_random(this, varargin)
            solver_opt = struct('method', 'halton',...
                'quasi_rand_seed', 1,...
                'num_quasi_rand_samples', 100 ...   % arbitrary - should be dim-dependant?  
            );
            solver_opt= varargin2struct_breach(solver_opt, varargin{:});

        
            this.solver = 'quasi_random';
            this.solver_options = solver_opt; 
        end

        function solver_opt = setup_random(this, varargin)
            solver_opt = struct( ...            
                'rand_seed', 1,...
                'num_rand_samples', 100 ...   % arbitrary - should be dim-dependant?  
            );
            solver_opt= varargin2struct_breach(solver_opt, varargin{:});
        
            this.solver = 'random';
            this.solver_options = solver_opt; 
        end

        
        function solver_opt = setup_corners(this, varargin)
            solver_opt = struct('num_corners', 100, ...   % arbitrary - should  be dim-dependant?  
            'group_by_inputs', true...
            );
            solver_opt= varargin2struct_breach(solver_opt, varargin{:});        
            
            this.solver = 'corners';
            this.solver_options = solver_opt; 
        end
 
        function solver_opt = setup_optimtool(this,varargin)
            solver_opt = optimset('Display', 'iter');
            this.display = 'off';
            solver_opt.lb = this.lb;
            solver_opt.ub = this.ub;
            this.solver = 'optimtool';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_fmincon(this,varargin)
            if this.verbose>=1
                disp('Setting options for fmincon solver');
            end
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
        
        function solver_opt = setup_fminsearch(this,varargin)
            %disp('Setting options for fminsearch solver');
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
        
        function solver_opt = setup_simulannealbnd(this,varargin)
            %disp('Setting options for simulannealbnd solver');
            solver_opt = saoptimset('Display', 'off');
            if this.max_time < inf
                solver_opt = saoptimset(solver_opt, 'MaxTime', this.max_time);
            end
            
            this.solver = 'simulannealbnd';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_cmaes(this,varargin)
            %disp('Setting options for cmaes solver - use help cmaes for details');
            solver_opt = cmaes();
            solver_opt.Seed = 0;
            solver_opt.LBounds = this.lb;
            solver_opt.UBounds = this.ub;
            solver_opt.start_sample = [];
            solver_opt.start_function_values = [];
            solver_opt= varargin2struct_breach(solver_opt, varargin{:});

            if isempty(this.x0)
               this.x0 = (this.ub+this.lb)/2; 
            end
            
            solver_opt.DispFinal='off';
            solver_opt.SaveVariables = 'off';
            solver_opt.LogModulo = 0;
            solver_opt.DispModulo = 0; % need to disable when running
            %multiple cmaes instances
            if this.use_parallel 
                solver_opt.EvalParallel = 'yes';
            end
            this.solver = 'cmaes';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_simulated_annealing(this, varargin)
            this.solver = 'simulated_annealing';
            solver_opt.lb = this.lb;
            solver_opt.ub = this.ub;
            solver_opt.start_sample = [];
            solver_opt.start_function_values = [];
            solver_opt= varargin2struct_breach(solver_opt, varargin{:});
            this.display = 'off';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_tomlab_glbfast(this, varargin)
            this.solver = 'tomlab_glbfast';
            solver_opt.lb = this.lb;
            solver_opt.ub = this.ub;
            this.display = 'off';
        end
        
        function solver_opt = setup_tomlab_lgo(this,varargin)
            this.solver = 'tomlab_lgo';
            solver_opt.lb = this.lb;
            solver_opt.ub = this.ub;
            this.display = 'off';
        end
        
        function solver_opt = setup_snobfit(this, varargin)
            % TODO: Put history here instead
            % (startSample, startFunctionValues)
            this.solver = 'snobfit';
            solver_opt.lb = this.lb;
            solver_opt.ub = this.ub;
            solver_opt.start_sample = [];
            solver_opt.start_function_values = [];
            solver_opt= varargin2struct_breach(solver_opt, varargin{:});   
            
            this.display = 'off';
            this.solver_options = solver_opt;
        end
		
		function solver_opt = setup_adaptive_corners(this, varargin)
            % Adaptive corners will shut down if no new falsification happens in
			% (relative_threshold)*(num_corners) consecutive simulations. 
            this.solver = 'adaptive_corners';
            solver_opt.lb = this.lb;
            solver_opt.ub = this.ub;
            solver_opt.relative_threshold = 0.25;
            solver_opt.num_corners = 100; % Arbitrary, make sure to change
            solver_opt= varargin2struct_breach(solver_opt, varargin{:});   
            
            this.display = 'off';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_uniform_random(this, varargin)
            this.solver = 'uniform_random';
            solver_opt.lb = this.lb;
            solver_opt.ub = this.ub;
            solver_opt.start_sample = [];
            solver_opt= varargin2struct_breach(solver_opt, varargin{:}); 
            this.display = 'off';
            this.solver_options = solver_opt;
        end
        
        function solver_opt = setup_gnmLausen(this)
            %disp('Setting options for GNM Lausen solver - use help gnm for details');
%             solver_opt = saoptimset('Display', 'off');
%             solver_opt.Seed = 0;
%             solver_opt.LBounds = this.lb;
%             solver_opt.UBounds = this.ub;
            
            solver_opt.maxRestarts = 15; % maximum (probablistic or degenerated) restarts
            solver_opt.maxEvals = 2500;  % maximum function evaluations
            solver_opt.nPoints = 5;      % number of random points per restart

            solver_opt.maxIter=250;      % maximum iterations per restart
            solver_opt.alpha = 1;        % reflection coeff
            solver_opt.beta = 0.5;       % contraction coeff
            solver_opt.gamma = 2;        % expansion coeff
            solver_opt.epsilon = 1e-9;   % T2 convergence test coefficient
            solver_opt.ssigma = 5e-4;    % small simplex convergence test coefficient
        
            if isempty(this.x0)
               this.x0 = (this.ub+this.lb)/2;
            end
            if ~isrow(this.x0)
                 solver_opt.xinit = this.x0';
            else
                solver_opt.xinit = this.x0;
            end
            
            this.solver = 'gnmLausen';
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
                
                case 'random'
                    res = this.solve_random();

                case 'quasi_random'
                    res = this.solve_quasi_random();
    
                case 'corners'
                    res = this.solve_corners();
                    
                case 'morris'
                    res = this.solve_morris();
                    
                case 'nelder_mead'
                    res = this.solve_nelder_mead();
                    
                case 'global_nelder_mead'
                    res = this.solve_global_nelder_mead();
                    
                case 'simulated_annealing'
                    % Simulated Annealing from S-TaLiRo
                    
                    
                    inputRanges = [this.lb this.ub];
                    
                    fun = this.objective;
                    
                    startSample = this.solver_options.start_sample;
                    
                    if isempty(startSample)
                        % No startSample given
                        [res, ~, startSample] = testron_SA(inputRanges, fun, this);
                    else
                        [res, ~, startSample] = testron_SA(inputRanges, fun, this, startSample);
                    end
                    
                case 'tomlab_glbfast'
                    % glbFast from TOMLAB
                    % Requires TOMLAB to be installed on the computer, plus
                    % a valid license file in the TOMLAB directory
                    
                    % To find example how to setup solver, see
                    % tomlab/quickguide/glbQG.m
                    if nargin < 2
                        % No startSample given
                        startSample = testronGetNewSample([this.lb this.ub]);
                    end
                    Name  = 'phi';
                    x_L   = this.lb;  % Lower bounds for x.
                    x_U   = this.ub;  % Upper bounds for x.
                    x_opt = [];
                    x_0   = startSample;
                    f_opt = [];
                    f_Low = [];             % Lower bound on function.
                    x_min = []; % For plotting
                    x_max = []; % For plotting

                    
                    fun = @(x, ~) this.objective(x);
                    
                    Prob  = glcAssign('tomlab_wrapper', x_L, x_U, Name, [], [], [], ...
                        [], [], [], x_0, ...
                        [], [], [], [], ...
                        f_Low, x_min, x_max, f_opt, x_opt);
                    
                    %Prob.PriLevOpt = 0; % No printing
                    Prob.optParam.MaxFunc = this.max_obj_eval;
                    Prob.optParam.fGoal = -eps;
                    
                    Prob.brProblem = this;
                    
                    Result = tomRun('glbFast', Prob, 0);
                    res = struct('bestRob',[],'bestSample',[],'nTests',[],'bestCost',[],'paramVal',[],'falsified',[],'time',[]);
                    res.bestSample = Result.x_k(:,1);
                    res.bestRob = Result.f_k;
                    
                case 'tomlab_lgo'
                    % LGO from TOMLAB
                    % Requires TOMLAB to be installed on the computer, plus
                    % a valid license file in the TOMLAB directory
                    
                    % To find example how to setup solver, see
                    % tomlab/quickguide/glbQG.m
                    if nargin < 2
                        % No startSample given
                        startSample = testronGetNewSample([this.lb this.ub]);
                    end
                    Name  = 'phi';
                    x_L   = this.lb;  % Lower bounds for x.
                    x_U   = this.ub;  % Upper bounds for x.
                    x_opt = [];
                    x_0   = startSample;
                    f_opt = [];
                    f_Low = [];             % Lower bound on function.
                    x_min = []; % For plotting
                    x_max = []; % For plotting
                    
                    fun = @(x, ~) this.objective(x);
                    
                    Prob  = glcAssign('tomlab_wrapper', x_L, x_U, Name, [], [], [], ...
                        [], [], [], x_0, ...
                        [], [], [], [], ...
                        f_Low, x_min, x_max, f_opt, x_opt);
                    
                    %Prob.PriLevOpt = 0; % No printing
                    Prob.optParam.MaxFunc = this.max_obj_eval;
                    Prob.optParam.fGoal = -eps;
                    
                    Prob.brProblem = this;
                    
                    Result = tomRun('lgo', Prob, 0);
                    res = struct('bestRob',[],'bestSample',[],'nTests',[],'bestCost',[],'paramVal',[],'falsified',[],'time',[]);
                    res.bestSample = Result.x_k(:,1);
                    res.bestRob = Result.f_k;
                    
                case 'snobfit'
                    res = this.solve_snobfit();
                    
                case 'cmaes'
                    %% set default sigma -- thao
                    if isempty(this.insigma)
                       sigmadefault = (this.ub-this.lb)/3;
                       if ~isrow(sigmadefault)
                           sigmadefault = sigmadefault'; 
                       end
                       this.insigma = sigmadefault';
                    end
                    %%
                    [x, fval, counteval, stopflag, out, bestever] = cmaes(this.objective, this.x0', this.insigma, this.solver_options);
                    res = struct('x',x, 'fval',fval, 'counteval', counteval,  'stopflag', stopflag, 'out', out, 'bestever', bestever);
                    this.add_res(res);
                    
                 case 'gnmLausen'
                    [x, fval, output, used_options] = gbnm(this.objective,this.lb,this.ub,this.solver_options);     
                    res = struct('x', output.usedPoints, 'fval', output.usedVals,...
                        'counteval', output.nEval, 'output', output, 'used_options', used_options);
%                     === example of output for 2D problem:
%                     usedPoints: [2×30 double]
%                     usedVals: [1×30 double]
%                     usedSimplex: {1×30 cell}
%                     reason: {1×15 cell}
%                     nEval: 831
                    this.add_res(res);
                    

                case 'meta'
                    res = this.solve_meta();        
                    
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
                    this.add_res(res);                    
            
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
                    this.add_res(res);
    
                case 'simulannealbnd'
                    % init seed if needed
                    if isfield(this.solver_options,'rand_seed')
                        rng(this.solver_options.rand_seed);
                    end
                    
                    this.display_status_header();
                    if this.use_parallel
                        warning('Parallel Computation not yet supported for Simulated Annealing.');
                        [x,fval,exitflag,output] = feval(this.solver, problem);
                        res = struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);                        
                        
%                         num_works = this.BrSys.Sys.Parallel;
%                         for idx = 1:num_works
%                             F(idx) = parfeval(this.solver, 4, problem);
%                         end
%                         res = cell(1,num_works);
%                         for idx = 1:num_works
%                             [completedIdx, x, fval,exitflag,output] = fetchNext(F);
%                             this.nb_obj_eval = this.nb_obj_eval + output.funccount;
%                             res{completedIdx} = struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);
%                         end
                    else
                        [x,fval,exitflag,output] = feval(this.solver, problem);
                        res = struct('x',x,'fval',fval, 'exitflag', exitflag, 'output', output);                        
                    end
                    this.add_res(res);
                case 'optimtool'
                    problem.solver = 'fmincon';
                    optimtool(problem);
                    res = [];
                    return;

                case 'binsearch'
                    res = solve_binsearch(this);
                    this.add_res(res);
                    
                case 'uniform_random'
                    % Initialize variables used
                    xbest = [];
                    fbest = Inf;
                    
                    printFlag = 1;
                    
                    for iterationCounter = 1:this.max_obj_eval
                        % Sample from uniform random distribution between
                        % this.lb and this.ub
                        x = (this.ub - this.lb).*rand(size(this.lb)) + this.lb;
                        
                        % Calculate robustness
                        rob = this.objective(x);
                        
                        if iterationCounter == 1
                            if numel(rob) == 1
                                disp(['Initial robustness value: ' num2str(rob)]);
                            else
                                disp(['Initial robustness value (using MAX semantics: ' num2str(min(rob))]);
                            end
                        end
                        
                        
                        
                        % Store if it's best
                        if min(rob) < fbest
                            xbest = x;
                            fbest = min(rob);
                            
                            if iterationCounter > 1
                                disp([num2str(iterationCounter) ': NEW BEST: ' num2str(fbest)]);
                            end
                        end
                        
                        if mod(iterationCounter, 10)==0 && ~this.stopping
                            fprintf([num2str(iterationCounter) ': Rob: ' num2str(min(rob)) '\t\tBEST:' num2str(fbest) '\n']);
                        end
                        
                        % Exit if robustness negative
                        if this.StopAtFalse && min(rob) < 0 
                            disp(['FALSIFIED at sample ' num2str(iterationCounter) '!']);
                            break
                        end
                    end
                    res = struct('bestRob',[],'bestSample',[],'nTests',[],'bestCost',[],'paramVal',[],'falsified',[],'time',[]);
                    res.bestSample = xbest;
                    res.bestRob = fbest;
                    
                case 'adaptive_corners'
                    res = this.solve_adaptive_corners();
                    
                otherwise
                    res = feval(this.solver, problem);
                    this.add_res(res);
                    
            end
            % TESTRON: Mute these outputs
            this.DispResultMsg(); 
            %this.Display_Best_Results(this.obj_best, this.x_best);
            
            %% Saving run in cache folder
            this.SaveInCache();
            
            % For some solvers we do not return the startSample
            if ~exist('startSample', 'var')
                startSample = '';
            end
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
        
        function X0 = init_basic_X0(this, n_samples)
            % returns initial vectors
            BrQ = this.BrSet.copy();
            BrQ.ResetParamSet();
            BrQ.SetParamRanges(this.params, [this.lb this.ub])
            BrC = BrQ.copy();
            nb_corners = this.solver_options.nb_max_corners;
            qseed = this.solver_options.quasi_rand_seed;
            nb_samples = this.solver_options.nb_new_trials;
            step = this.solver_options.start_at_trial;
            
            BrC.P = CreateParamSet(BrC.Sys,this.params,[this.lb this.ub]);
            % JOHAN CHANGE
            if numel(this.params) < 5
                % Standard case, use CornerSample
                fprintf('%d varying parameters, using standard CornerSample\n',numel(this.params));
                BrC.CornerSample();
            else
                % Too many parameter combinations to enumerate
                % Use QuasiRandomSample instead
                if ~exist('n_samples', 'var')
                    n_samples = 180; % Arbitrary choice
                end
                fprintf('%d varying parameters, using QuasiRandomSample(%d) (TestronRefine) to not run out of memory\n',numel(this.params),n_samples);
                BrC.QuasiRandomSample(n_samples);
                
                % Change nb_samples, which is a "Breach" variable
                nb_samples = n_samples;
            end
            % END JOHAN CHANGE
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
                    BrQ.QuasiRandomSample(qnb_samples, qseed);
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
            
            if numel(this.Spec.req_monitors)>1
                fun_obj = @(x)(min(this.objective(x),[],1)); % for basic multi-objective support                
            else
                fun_obj = this.objective;
            end
            
            
            problem = struct('objective', fun_obj, ...
                'fitnessfcn', fun_obj, ... % for ga
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
        
        % check the MIP options for supported solvers and change the this.params
        function setup_mixed_int_optim(this, method)
            if ~exist('method')||isempty(method)
                method = 'map_enum_to_int';
            end
            
            switch method
                case 'map_enum_to_int'
                    if ismember(this.solver, this.mixed_integer_optim_solvers)
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
                    end
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
        
        %% Objective function and wrapper        
        
        function [obj, cval, x_stoch] = objective_fn(this,x)
            
            % reset this.Spec
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
            
            % For falsification, default objective_fn is mostly robust satisfaction of the least
            this.robust_fn(x);
            robs = this.Spec.traces_vals;

            if (~isempty(this.Spec.traces_vals_precond))
                num_tr = size(this.Spec.traces_vals_precond,1);
                precond_robs = zeros(num_tr,1);
                for itr = 1:num_tr
                    precond_robs(itr) = min(this.Spec.traces_vals_precond(itr,:));
                    if  precond_robs(itr)<0
                        robs(itr,:)= -precond_robs(itr);
                    end
                end
            end
            NaN_idx = isnan(robs); % if rob is undefined, make it inf to ignore it
            robs(NaN_idx) = inf;
            obj = min(robs,[],1)';
            cval = inf;
            if (~isempty(this.Spec.traces_vals_precond))
                NaN_idx = isnan(precond_robs); % if rob is undefined, make it inf to ignore it
                precond_robs(NaN_idx) = inf;
                cval = min(precond_robs,[],1)';
            end
            
                        
        end
        
        function [fval, cval, x_stoch] = objective_wrapper(this,x)
            % reset this.Spec
            % objective_wrapper calls the objective function and wraps some bookkeeping
            if size(x,1) ~= numel(this.params)
                x = x';
            end
            
            nb_eval =  size(x,2);
            fval = inf*ones(size(this.Spec.req_monitors,2), nb_eval);            
            cval = inf*ones(1, nb_eval);
            x_stoch = nan(numel(this.stochastic_params), nb_eval);
            
            fun = @(isample) this.objective_fn(x(:, isample));
            nb_iter = min(nb_eval, this.max_obj_eval);
     
            if this.stopping()==false
                if nb_iter == 1 || ~this.use_parallel
                    iter = 0;
                    while ~this.stopping()&&iter<size(x,2)
                        iter= iter+1;
                        
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
                            [fval(:,iter), cval(:,iter), x_stoch(:,iter)] = fun(iter);
                            
                            % Normalize the function value based on average
                            % robustness stated
                            fval(:, iter) = fval(:, iter)./this.avgRobForNormalization;
                            
                            % logging and updating best
                            this.time_spent = toc(this.time_start);
                            this.LogX(x(:, iter), fval(:,iter), cval(:,iter), x_stoch(:,iter));

                            % update status
                            if ~rem(this.nb_obj_eval,this.freq_update)
                                this.display_status();
                                % callback
                                if ~isempty(this.callback_obj)
                                   e.name ='obj_computed';
                                   e.values.fval = fval;                                   
                                   e.values.nb_obj_eval=this.nb_obj_eval;
                                   this.callback_obj(this, e);
                                end
                            end

                        end
                                                                        
                    end
                else % Parallel case
                                        
                    for batch_counter = 1:this.parallelBatchSize:nb_iter                                                    
                        start_index = batch_counter;
                        end_index = min(batch_counter + this.parallelBatchSize - 1, nb_iter);
                                                
                        % Launch tasks                        
                        for iter = start_index:end_index                            
                            par_f(:,iter) = parfeval(fun,3, iter);                            
                        end
                        
                        num_this_batch = end_index-start_index+1;                                                
                        for idx = 1:num_this_batch
                            
                            [completedIdx(idx), fval_batch, cval_batch ,x_stoch_batch] = fetchNext(par_f);
                            
                            % Timing to get new fetch
                            this.time_spent = toc(this.time_start);
                            
                            % Normalize the function value based on average
                            % robustness stated
                            fval_batch = fval_batch./this.avgRobForNormalization;                                                        
                         
                            % we log in the order of arrival, but will 
                            % reorder after the batch using completedIdx
                            this.LogX(x(:,completedIdx(idx)),fval_batch,cval_batch,x_stoch_batch); 
                            
                            % this way we can update status
                            if ~rem(this.nb_obj_eval,this.freq_update)
                                this.display_status(fval_batch, cval_batch);
                            end
                            
                            fval(:,completedIdx(idx)) = fval_batch;
                            cval(:,completedIdx(idx)) = cval_batch;
                            x_stoch(:,completedIdx(idx)) = x_stoch_batch;
                            
                            if this.stopping()
                                cancel(par_f);
                                break
                            end
                            
                        end
                            
                        % Let's reorder the log 
                        [~, ia] = sort(completedIdx);
                        this.obj_log(:,end-numel(completedIdx)+1:end) = this.obj_log(:,ia);
                        this.X_log(:,end-numel(completedIdx)+1:end) = this.X_log(:,ia);
                        this.X_stochastic_log(:,end-numel(completedIdx)+1:end) = this.X_stochastic_log(:,ia);  
                                                                                                
                    end
                end
            else
                fval = this.obj_best*ones(1, nb_eval);
                
            end                       
            
        end
        
        function b = stopping(this)
            b = (this.time_spent >= this.max_time) ||...
                (this.nb_obj_eval>= this.max_obj_eval) || ...
                (this.num_consecutive_constraints_failed >= this.max_consecutive_constraints_failed)||...
                this.is_paused;
        end
              
        %% Misc methods
        function x = CheckinDomain(this,x)
          for ip = 1:numel(this.params)
                x(ip) = this.domains(ip).checkin(x(ip));  
          end
        end
        
        function LogX(this, x, fval, cval, x_stoch)
            % LogX logs values tried by the optimizer

            if cval>=0
                this.num_consecutive_constraints_failed = 0;
                
                this.X_log = [this.X_log x];
                this.obj_log = [this.obj_log fval];
                this.X_stochastic_log = [this.X_stochastic_log x_stoch];
                if (this.log_traces)%&&~(this.use_parallel)&&~(this.BrSet.UseDiskCaching) % FIXME - logging flags and methods need be revised
                    if isempty(this.BrSet_Logged)
                        this.BrSet_Logged = this.BrSys.copy();
                        this.R_log = this.Spec.copy();
                    else
                        this.BrSet_Logged.Concat(this.BrSys, true); % fast concat
                        this.R_log.Concat(this.Spec, true);
                    end
                end
                
                [fmin , imin] = min(min(fval));
                x_min =x(:, imin);
                if fmin < this.obj_best
                    this.x_best = x_min;
                    this.obj_best = fval(:,imin);
                    this.BrSet_Best = this.BrSys.copy(); % could be more efficient - also suspicious when several x...
                end
                
                % num_eval
                this.nb_obj_eval= numel(this.obj_log(1,:));
            else
                 this.num_consecutive_constraints_failed = this.num_consecutive_constraints_failed+1;
                 this.num_constraints_failed = this.num_constraints_failed+1;                
            end
            if rem(this.nb_obj_eval+this.num_constraints_failed,this.freq_update)
                % TODO: what should be here?
            end
        end
         
        function [BrOut, Berr, BbadU] = GetBrSet_Logged(this)
            % GetBrSet_Logged gets BreachSet object containing parameters and traces computed during optimization
            BrOut = this.BrSet_Logged;
            if isempty(BrOut)||(~BrOut.hasTraj())
                BrOut = this.BrSys.copy();
                BrOut.ResetSimulations();               
                if isempty(this.stochastic_params) 
                    BrOut.SetParam(this.params, this.X_log, 'combine');
                else
                    BrOut.SetParam([this.params this.stochastic_params], [this.X_log ; this.X_stochastic_log], 'combine');
                end
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
            if isempty(this.X_log)
                Rlog =[];
            else
                Rlog = this.GetBrSet_Logged(varargin{:});
            end
        end
        
        function Rbest = GetBest(this,varargin)
            Rbest = this.GetBrSet_Best(varargin{:});
        end
    
        
        function Display_X(this, param_values)
            if ~isempty(param_values)
                
                for ip = 1:numel(this.params)
                    value = param_values(ip);
                    fprintf( '        %s = %g\n', this.params{ip},value)
                end
                fprintf('\n');
            end
        end

        function DispResultMsg(this)
            if ~strcmp(this.display, 'off')
                % display status of last eval if not already done
                if rem(this.nb_obj_eval,this.freq_update)
                    this.display_status();
                end
                
                % DispResultMsg message displayed at the end of optimization
                if this.time_spent >= this.max_time
                    fprintf('\n Stopped after max_time was reached.\n');
                end
                
                if this.nb_obj_eval >= this.max_obj_eval
                    fprintf('\n Stopped after max_obj_eval was reached (maximum number of objective function evaluations).\n' );
                end
                
                if this.num_consecutive_constraints_failed >= this.max_consecutive_constraints_failed 
                    fprintf('\n Stopped after maximum number of consecutive constraints \n(max_consecutive_constraints_failed=%g) failed was reached.\n', this.max_consecutive_constraints_failed);                
                end
            end
        end
        
        function Display_Best_Results(this, best_fval, param_values)
            if ~strcmp(this.display, 'off')
                if ~isempty(param_values)
                    fprintf('\n ---- Best value %g found with\n', min(best_fval));
                    this.Display_X(param_values);
                else
                    fprintf('\n ---- Failed to find parameter values satisfying constraints.\n\n');
                end
            end
        end
   
    end
    
    methods (Access=protected)
        
        function display_status_header(this)
            if ~isempty(this.Spec.precond_monitors)
                hd_st = sprintf(  '#calls (max:%5d)        time spent (max: %g)     [current  obj]     (current best)   [constraint (num failed)]\n',...
                    this.max_obj_eval, this.max_time);
            else
                hd_st = sprintf(  '#calls (max:%5d)        time spent (max: %g)     [current  obj]     (current best) \n',...
                    this.max_obj_eval, this.max_time);
            end
            fprintf(hd_st);
        end
        
        function display_status(this,fval,cval)
            
            if ~strcmp(this.display,'off')
                if nargin==1
                    fval = this.Spec.val; 
                    if ~isempty(this.Spec.precond_monitors)
                        cval = min(min(this.Spec.traces_vals_precond)); % bof bof
                    end
                end
                
                st__= sprintf('     %5d                    %7.1f               [%s]     (%s)', ...
                    this.nb_obj_eval, this.time_spent,  num2str(fval','%+5.5e '), num2str(this.obj_best', '%+5.5e '));
                if ~isempty(this.Spec.precond_monitors)
                    st__ = sprintf([st__ '           [%s (%g)]\n'], num2str(cval', '%+5.5e '), this.num_consecutive_constraints_failed);
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
            BrOut = []; 
            if isempty(B)
                BrOut = [];
                return;
            end
             B.Sys.Verbose=1;
           
            [idx_ok, idx_sim_error, idx_invalid_input, st_status]  = B.GetTraceStatus();
            
            if ~isempty(idx_sim_error)||~isempty(idx_invalid_input)
                [B, Berr, BbadU] = FilterTraceStatus(B);
                this.disp_msg(['Note: ' st_status],1);
            end
            
            if ~isempty(idx_ok)&&B.hasTraj()
                if isa(this.R_log, 'BreachRequirement')&&this.R_log.hasTraj()&&...
                        numel(idx_ok)==numel(this.R_log.P.traj)
                    BrOut=this.R_log;
                    BrOut.BrSet= B;
                else
                    BrOut = this.R0.copy();
                    BrOut.Eval(B);
                end            
            end
        end
 
        function add_res(this,res)
        % appends new optimization result
            if isempty(this.res)
                this.res = res;
            elseif isstruct(this.res)
                res_list = {this.res};
                this.res = res_list;
            else
                this.res{end+1} = res; 
            end
        end
    end
    
    
    
end
