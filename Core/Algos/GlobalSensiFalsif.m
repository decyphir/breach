classdef GlobalSensiFalsif < FalsificationProblem
    
    methods
        
        function this = GlobalSensiFalsif(BrSys, phi, params, ranges)
           
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
            this = this@FalsificationProblem(super_args{:});
            this.obj_best=inf;
        end
                
        function solver_opt = setup_morris(this, varargin)
            % solve_corners works with quasi-random sampling
            %
            
            solver_opt = struct('num_path', 100, ...   % number of paths, i.e., set of
                ...                    % samples providing 1 pair of samples per dim
                'size_grid', 5 ...     % number of grid levels, intervals in each dim
                );
            solver_opt = varargin2struct_breach(solver_opt, varargin{:});
            
            
            this.solver = 'morris';
            this.solver_options = solver_opt;
            
        end
        
        
        function res = solve(this)
            % Morris method - compute simulations and display resulting
            % sensitivities
            %
            
            r = this.solver_options.num_path;
            p = this.solver_options.size_grid;
            dim = numel(this.params);
            
            % number of samples is going to be (dim+1)*r
            num_samples = (dim+1)*r;
            if this.max_obj_eval<num_samples % we have to reevaluate the number of path
                r = floor(this.max_obj_eval/(dim + 1));
            end
            
            Sys= CreateSystem({},this.params, this.x0);
            P  = CreateParamSet(Sys, this.params, [this.lb, this.ub]);
            Pr = pRefine(P, p,r);
            X0 = Pr.pts;
            
            if ~strcmp(this.display,'off')
                fprintf('\n Running %g samples using Morris'' method\n', size(X0, 2));
                this.display_status_header();
            end
            
            res = this.FevalInit(X0);
            
            [res.mu, res.mustar, res.sigma, res.sigmastar, res.EE] = EEffects(res.fval, Pr.D, p);
            res.params = this.params;
            % this.add_res(res);
            
            if ~strcmp(this.display, 'off')
                this.display_morris_result(res);
            end
        end
        
        function display_morris_result(this, res, varargin)
            
            opt.SortBy = 'mu';
            opt = varargin2struct(opt,varargin{:});
            
            switch opt.SortBy
                case 'mu'
                    [~ , order] = sort(res.mu);
                case 'mustar'
                    [~ , order] = sort(res.mustar);
                case 'sigma'
                    [~ , order] = sort(res.sigma);
                case 'sigmastar'
                    [~ , order] = sort(res.sigmastar);
            end
            
            fprintf('mu        mustar    sigma     sigmastar\n')
            
            for ip = 1:numel(res.params)
                idx = order(ip);
                %format = '%+5.5e';
                format = '%+5.5G';
                
                st_format = [format '   ' format '   ' format '   ' format '   ---> %s\n' ];
                fprintf(st_format, ...
                    res.mu(idx), res.mustar(idx), res.sigma(idx), res.sigmastar(idx), res.params{idx})
            end
            
            
            
            
            
        end
        
        
    end
    
    
end

