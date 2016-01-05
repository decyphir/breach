classdef BreachSystem < BreachSet
    %BreachSystem    Defines a simplified API for Breach functionalities.
    %                It combines a system structure (Sys) and parameter set (P)
    %                into one object so that basic operations can be done in one
    %                command instead of several and with fewer arguments.
    
    properties
        Sys       % Legacy Breach system structure
        Specs     % A set (map) of STL formulas
        ParamRanges
        SignalRanges
        Problem
    end
    
    methods
        
        %% Constructor
        function this = BreachSystem(varargin)
            this.Specs = containers.Map();
            
            switch nargin
                case 0 % do nothing
                case 1 % Should be a Sys structure
                    
                    inSys = varargin{1};
                    if isaSys(inSys)
                        this.Sys = inSys;
                        this.ResetParamSet();
                    else
                        error('BreachObject with one argument assumes that the argument is a system structure.')
                    end
                otherwise % creates an extern system
                    if ~(exist(varargin{1})==4) % tests if the first argument is a Simulink model
                        this.Sys = CreateExternSystem(varargin{:});
                    else
                        warning('First argument is the name of a Simulink model - consider using BreachSimulinkSystem instead');
                    end
            end
            
            if isaSys(this.Sys) % Note: we ignore initial conditions for now in ParamRanges
                % OK for Simulink, less so for ODEs...
                this.ParamRanges = [this.Sys.p(this.Sys.DimX+1:end) this.Sys.p(this.Sys.DimX+1:end)];
                this.SignalRanges = [];
                this.ResetParamSet();                        
            end            
        end
        
        %% Parameters
        % ResetParamSet Reset parameter set based on ParamRanges
        function ResetParamSet(this)
            ipr = find(diff(this.ParamRanges'));
            ranges =this.ParamRanges(ipr,:);
            if (isempty(ipr))
                this.P = CreateParamSet(this.Sys);
                this.P.epsi(:,:)=0;
            else
                this.P = CreateParamSet(this.Sys, this.P.ParamList(ipr+this.Sys.DimX),ranges);
            end
            
        end
        
        % Get and set default parameter values (defined in Sys)
        function values = GetDefaultParam(this, params)
            values = GetParam(this.Sys,params);
        end
        function SetDefaultParam(this, params, values)
            this.Sys = SetParam(this.Sys,params, values);
        end
               
        %% Signals plots and stuff        
        function SetTime(this,tspan)
           this.Sys.tspan = tspan; 
        end
        function time = GetTime(this,tspan)
           time = this.Sys.tspan;
        end
         
        % Performs a simulation from the parameter vector(s) defined in P
        function Sim(this,tspan)
            if nargin==1
                tspan = this.Sys.tspan;
            end
            this.P = ComputeTraj(this.Sys, this.P, tspan);
        end
        
        % Update ranges for variables from trajectories in P
        function val = UpdateSignalRanges(this)
            
            if isfield(this.P, 'traj')
                if isempty(this.SignalRanges)
                    this.SignalRanges = ones(this.Sys.DimX,2);
                    minX = +inf*ones(this.Sys.DimX,1);
                    maxX = -inf*ones(this.Sys.DimX,1);
                else
                    minX = this.SignalRanges(:,1);
                    maxX = this.SignalRanges(:,2);
                end
                val=inf;
                for itraj = 1:numel(this.P.traj)
                    traj = this.P.traj(itraj);
                    traj_maxX = max(traj.X,[], 2);
                    traj_minX = min(traj.X,[], 2);
                    dist_maxX = min(maxX-traj_maxX);
                    dist_minX = min(traj_minX-minX);
                    val= min( [val dist_maxX dist_minX] );
                    minX = min([traj_minX minX],[],2);
                    maxX = max([traj_maxX maxX],[],2);
                end
                this.SignalRanges = [minX, maxX];
            end
        end
        
        % Get signal names
        function SigNames = GetSignalNames(this)
            SigNames = this.Sys.ParamList(1:this.Sys.DimX);
        end
        
        % Plot signals
        function h = PlotSignals(this, varargin)
            Sim(this);
            figure;
            h = SplotVar(this.P, varargin{:});
        end
        
        %% Specs
        % Add (a) specs
        function phi = AddSpec(this, varargin)
            global BreachGlobOpt
            if nargin==2
                if strcmp(class(varargin{1}),'STL_Formula');
                    phi = varargin{1};
                elseif ischar(varargin{1})
                    phi_id = MakeUniqueID([this.Sys.name '_spec'],  BreachGlobOpt.STLDB.keys);
                    phi = STL_Formula(phi_id, varargin{:});
                end
            end
            
            % checks signal compatibility
            [~,sig]= STL_ExtractPredicates(phi);
            i_sig = FindParam(this.Sys, sig);
            sig_not_found = find(i_sig>this.Sys.DimP, 1);
            if ~isempty(sig_not_found)
                error('Some signals in specification are not part of the system.')
            end
            
            this.Specs(get_id(phi)) = phi;
        end
        
        function val = CheckSpec(this, spec)
            if ~exist('spec','var')
                spec = this.Specs.values
            else
                if iscell(spec)
                    for cur_spec = spec
                        this.AddSpec(cur_spec{1});
                    end
                else
                    spec = this.AddSpec(spec);
                end
                
                [this.P, val] = SEvalProp(this.Sys,this.P,spec);
            end
        end
        
        % Monitor spec on reference traj
        function [rob, tau] = GetRobustSat(this, phi, params, values)
            if nargin==1
                phi = this.spec;
                params = {};
                values = [];
            end
            
            if nargin==2
                params = {};
                values = [];
            end
            
            if ischar(phi)
                this__phi__ = STL_Formula('this__phi__', phi);
            else
                this__phi__ = phi;
            end
            
            if ~isempty(params)
                this.P = SetParam(this.P, params, values);
            end
            
            Sim(this);
            [rob, tau] = STL_Eval(this.Sys, this__phi__, this.P, this.P.traj,0);
        end
        
        % Return a function of the form robfn: p -> rob such that p is a
        % vector of values for parameters and robfn(p) is the
        % corresponding robust satisfaction
        function robfn = GetRobustSatFn(this, phi, params)
            
            if ischar(phi)
                this__phi__ = STL_Formula('this__phi__', phi);
                robfn = @(values) GetRobustSat(this, this__phi__, params, values);
            else
                robfn = @(values) GetRobustSat(this, phi, params, values);
            end
            
        end
        
        %% Falsification of properties
        function FalsifySpec(this, phi,falsif_opt, param_prop)
            % See Falsify
            this.P = Falsify(this.Sys, phi, falsif_opt, param_prop);
            
        end
        
        %% Sensitivity analysis
        % FIXME interface not complete
        function SensiSpec(this,phi)
            this.ResetParamSet();
            opt.tspan = 0:0.1:50;
            opt.params = this.P.dim;
            opt.lbound = this.ParamRanges(this.P.dim-this.P.DimX,1)';
            opt.ubound = this.ParamRanges(this.P.dim-this.P.DimX,2)';
            opt.plot = 2;
            
            [mu, mustar, sigma, Pr]= SPropSensi(this.Sys, this.P, phi, opt);
            this.P=Pr;
        end
        
        %% Mining
        function [p, rob] = MineSpec(this, phi, falsif_opt, prop_opt, iter_max)
            [p, rob, Pr] = ReqMining(this.Sys, phi, falsif_opt, prop_opt, iter_max);
            this.P = Pr;
        end
        
        function [res]  = MaxSatSpec(this, phi, params, problem)
            
            if ~exist('params','var')
                params = this.P.ParamList(this.P.dim);
            end
            
            if ischar(params)
                params = {params};
            end
            
            % Robust satisfaction function to maximize
            robust_fn = this.GetRobustSatFn(phi, params);
            
            % Paramter ranges
            ranges = this.GetParamRanges(params);
            lb = ranges(:,1);
            ub = ranges(:,2);
            
            % if range is singular, assumes unconstrained
            issame  = find(ub-lb==0);
            lb(issame) = -inf;
            ub(issame) = inf;
            
            % Initial value
            x0 = this.GetParam(params)';
            if isempty(x0)
                default_params = get_params(phi);
                for ip = 1:numel(params)
                    x0(ip) = default_params.(params{ip});
                end
            end
            
            % Problem structure
            if ~exist('problem','var')
                problem = BreachProblem(robust_fn, x0, 'max');
            else
                problem.robust_fn = robust_fn;
                problem.objective = @(x) (problem.max_robust_fn(x));
                problem.x0 = x0;
            end
            
            problem.lb = lb;
            problem.ub = ub;
            res = problem.solve();
            this.Problem = problem;
        end
        
        % Generic param synthesis function
        function [xopt] = MineParam(this, phi, params, problem)
            % Syntax BrObj.MineParam(phi,phi, params, options)
            
            if ~exist('params','var')
                params = this.P.ParamList(this.P.dim)
            end
            
            if ischar(params)
                params = {params};
            end
            
            % Robust satisfaction function to maximize
            robust_fn = this.GetRobustSatFn(phi, params);
            
            % Paramter ranges
            ranges = this.GetParamRanges(params);
            lb = ranges(:,1);
            ub = ranges(:,2);
            
            % if range is singular, assumes unconstrained
            issame  = find(ub-lb==0);
            lb(issame) = -inf;
            ub(issame) = inf;
            
            % Initial value
            x0 = this.GetParam(params)';
            if isempty(x0)
                default_params = get_params(phi);
                for ip = 1:numel(params)
                    x0(ip) = default_params.(params{ip});
                end
            end
            
            % Problem structure
            if ~exist('problem','var')
                problem = BreachProblem(robust_fn, x0, 'tight');
            else
                problem.robust_fn = robust_fn;
                problem.objective = @(x) (problem.tight_robust_fn(x));
                problem.x0 = x0;
            end
            
            problem.lb = lb;
            problem.ub = ub;
            xopt = problem.solve();
            
        end
        
        %% Printing
        function PrintSignals(this)
            if isempty(this.SignalRanges)
                disp( 'Signals:')
                disp( '-------')
                for isig = 1:this.Sys.DimX
                    fprintf('%s\n', this.Sys.ParamList{isig});
                end
            else
                
                fprintf('Signals (in range estimated over %d simulations):\n', numel(this.P.traj))
                disp('-------')
                for isig = 1:this.Sys.DimX
                    fprintf('%s in  [%g, %g]\n', this.Sys.ParamList{isig}, this.SignalRanges(isig,1),this.SignalRanges(isig,2));
                end
            end
            disp(' ')
        end
        
        function PrintParams(this)
            nb_pts= this.GetNbParams();
            if (nb_pts==1)
                disp('Parameters:')
                disp('----------')
                for ip = this.Sys.DimX+1:numel(this.P.ParamList)
                    fprintf('%s=%g',this.P.ParamList{ip},this.P.pts(ip,1));
                    rg = this.ParamRanges(ip-this.Sys.DimX,2)-this.ParamRanges(ip-this.Sys.DimX,1);
                    if rg>0
                        fprintf(', can vary in [%g, %g]',this.ParamRanges(ip-this.Sys.DimX,1),this.ParamRanges(ip-this.Sys.DimX,2));
                    end
                    fprintf('\n');
                end
            else
                fprintf('Parameters (%d values):\n',nb_pts);
                disp('-------------------------');
                for ip = this.Sys.DimX+1:numel(this.P.ParamList)
                    rg = this.ParamRanges(ip-this.Sys.DimX,2)-this.ParamRanges(ip-this.Sys.DimX,1);
                    if rg>0
                        fprintf('%s',this.P.ParamList{ip});
                        fprintf(' varying in [%g, %g]\n',this.ParamRanges(ip-this.Sys.DimX,1),this.ParamRanges(ip-this.Sys.DimX,2));
                    else
                        fprintf('%s=%g\n',this.P.ParamList{ip},this.P.pts(ip,1));
                    end
                end
            end
            
            disp(' ')
        end
        
        function PrintSpecs(this)
            disp('Specifications:')
            disp('--------------')
            keys = this.Specs.keys;
            for is = 1:numel(keys)
                fprintf('%s\n',keys{is});
            end
            disp(' ');
        end
        
        function disp(this)
            disp(['BreachObject interfacing model ' this.Sys.name '.']);
        end
        
        %% GUI
        function RunGUI(this)
            %assignin('base','P__', this.P);
            %assignin('base','Sys__', this.Sys);
            
            %evalin('base', 'Breach(Sys__,''P__'');');
            P.P = this.P;
            phis=  this.Specs.keys;
            
            Psave(this.Sys, 'Pthis', this.P);
            
            if (~isempty(phis))
                Propsave(this.Sys, phis{:});
            end
            Breach(this.Sys);
        end
        
        %% Misc
        % FIXME Resets the system to nominal parameters
        function Reset(this)
            this.P = CreateParamSet(this.Sys);
        end
        
        % Removes computed trajectories
        function ResetSimulations(this)
            this.P = SPurge(this.P);
        end
        
        % Make a copy of a handle object - works because no property is
        % itself a handle object.
        function new = copy(this)
            % Instantiate new object of the same class.
            new = feval(class(this));
            
            % Copy all non-hidden properties.
            p = fieldnames(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end
        
    end
end
