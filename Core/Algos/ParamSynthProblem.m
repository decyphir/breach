classdef ParamSynthProblem < BreachProblem
    % ParamSynthProblem Parameter Synthesis variant of BreachProblem - we minimize abs value of
    % robustness while keeping it positive.
    %
    % The default solver is binsearch, which assumes monotonicity wrt
    % parameters. The algorithm tries to infer monotonicity - if this
    % fails, monotonicity can be specified manually with the field monotony.
    % E.g.,
    %   paramsynth_pb.solver_options.monotony = [-1 1];
    % 
    % means the objective function is decreasing wrt the first parameter and 
    % increasing wrt the second parameter. 
    %
    %
    % ParamSynthProblem Properties
    %    BrSet_Synth   -  a BreachSet with synthesized parameters
    %    X_Synth       -  Synthesized parameters
    %
    % ParamSynthProblem Methods
    %   GetBrSet_Synth - returns BrSet_Synth
    %
    % See also BreachProblem
    
    properties
        BrSet_Synth
        X_Synth
        tightness
    end
    
    methods
        
        % Constructor calls parent constructor
        function this = ParamSynthProblem(BrSys, phi, params, ranges)
            
            Br = BrSys.copy();
               
            switch nargin
                case 2
                    sys_var = Br.GetSysVariables();
                    if ~isempty(sys_var)
                        Br.ResetDomain(sys_var);
                    end
                    
                    super_args{1} = Br;
                    super_args{2} = phi;
                    req_params = Br.GetReqVariables();
                    super_args{3} = req_params;
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
            this.constraints_fn = @(x) tight_constraints(this,x);
            this.tightness = (this.ub-this.lb)/100;  % arbitrary
            if ~isfield(this.BrSet.P, 'traj') 
                this.BrSet.Sim(); % ensures that we have trajectory to infer parameter from 
            end
            this.setup_binsearch();
        end
        
        % default objective function for synthesis
        function obj = objective_fn(this,x)
            this.Spec = this.R0.copy(); 
            rob = this.robust_fn(x);
            rob = min(rob);
            obj = max([rob 0]) + 1000*max([-rob 0]);
            
        end
        
        % Computes the tightness constraints - doesn't account for
        % monotonicity yet
        function const = tight_constraints(this,x)
            nb_params = numel(this.params);
            const_vec = zeros(nb_params,1);
            rob = min(this.robust_fn(x));
            for ip = 1:nb_params
                delta = this.tightness(ip);
                xi = x(ip);
                % if x close to boundary, ok
                if (xi+delta > this.ub(ip))||(xi-delta<this.lb(ip))
                    const_vec(ip) = 1;
                else
                    xplus = x;
                    xplus(ip) = x(ip)+ delta;
                    xminus = x;
                    xminus(ip) = x(ip)- delta;
                    robplus = min(this.robust_fn(xplus));
                    robminus = min(this.robust_fn(xminus));
                    
                    dplus = - robplus*rob;
                    dminus = - robminus*rob;
                    
                    const_vec(ip) = max([dplus dminus]);  
                
                end
            end
            const = max(const_vec);
        end
        
        function ResetObjective(this, varargin)
            ResetObjective@BreachProblem(this, varargin{:});
            this.X_Synth = [];
            this.BrSet_Synth = [];
            this.x0 = this.x0(:,1); % need double checking... 
        end
        
        
        function BrSynth = GetBrSet_Synth(this)
            BrSynth = this.BrSet_Best.copy();
            BrSynth.P = Sselect(BrSynth.P,1); % can do better...
        end
        
        function [X_Synth, res] = solve(this)
            res = solve@BreachProblem(this);
            this.X_Synth = this.x_best;
            X_Synth = this.X_Synth;
        end
        
        function DispResultMsg(this)
            fprintf('\n ---- Best objective value found: %g, with\n', this.obj_best);
            param_values = this.x_best;
            for ip = 1:numel(this.params)
                fprintf( '        %s = %g\n', this.params{ip},param_values(ip))
            end
            fprintf('\n');
            
        end
        
    end
end