classdef BreachSystem < BreachSet
    %BreachSystem  A simplified API class for Breach.
    %
    % It combines a system structure (Sys) with a parameter set (P)
    % into one object so that basic operations can be done in one
    % command instead of several and with fewer arguments. BreachSystem
    % class is derivated from BreachSet. Type help BreachSet to view
    % properties and methods of the parent class.
    %
    % BreachSystem Properties
    %   Specs  - a set of Signal Temporal Logic (STL) formulas.
    %
    %
    % BreachSystem methods
    %   Sim           - Simulate the system for some time using every parameter vectors.
    %   SetTime       - Set a default time for simulation (array or end time)
    %   GetTime       - Get default time for simulation
    %   AddSpec       - Add a specification to the system
    %   CheckSpec     - Checks (return robust satisfcation of) a given specification or all added specs
    %   PlotRobustSat - Plots robust satisfaction of a given STL formula against time
    %   PlotRobustMap - Plots (1d or 2d) robust satisfaction against parameter values
    %   RunGUI        - Open Breach legacy GUI, allowing for interactively exploring parameters, traces and specifications
    %
    %See also BreachSet.
    
    properties
        Sys       % Legacy Breach system structure
        Specs     % A set (map) of STL formulas
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
        % Get and set default parameter values (defined in Sys)
        function values = GetDefaultParam(this, params)
            values = GetParam(this.Sys,params);
        end
        function SetDefaultParam(this, params, values)
            this.Sys = SetParam(this.Sys,params, values);
        end
        function SetP(this,P)
            
            if isaP(P)
                this.P =P;
                this.UpdateParamRanges();
            else
                error('Argument should a Breach parameter structure.');
            end
            
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
        
        
        %% Specs
        % Add (a) specs % TODO double check
        function phi = AddSpec(this, varargin)
            global BreachGlobOpt
            if nargin==2
                if strcmp(class(varargin{1}),'STL_Formula');
                    phi = varargin{1};
                elseif ischar(varargin{1})
                    phi_id = MakeUniqueID([this.Sys.name '_spec'],  BreachGlobOpt.STLDB.keys);
                    phi = STL_Formula(phi_id, varargin{2});   % Can't be right ...
                end
            else
                return;
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
            end
            [this.P, val] = SEvalProp(this.Sys,this.P,spec);
        end
        
        % Monitor spec on reference traj
        function [rob, tau] = GetRobustSat(this, phi, params, values, tau)
            
            if nargin < 5
                tau = 0;
            end
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
        function [robfn, BrSys] = GetRobustSatFn(this, phi, params)
            
            BrSys = this.copy();
            
            if ischar(phi)
                this__phi__ = STL_Formula('this__phi__', phi);
                robfn = @(values) GetRobustSat(BrSys, this__phi__, params, values);
            else
                
                robfn = @(values) GetRobustSat(BrSys, phi, params, values);
            end
            
        end
        
        
        % Plots satisfaction signal
        function PlotRobustSat(this, phi, depth, tau, ipts)
            % check arguments
            if(~exist('ipts','var')||isempty(ipts))
                ipts = 1:size(this.P.pts,2);
            end
            
            if(~exist('tau','var')||isempty(tau)) % also manage empty cell
                tau = [];
            end
            
            if ~exist('depth','var')
                depth = inf;
            end
            
            figure;
            SplotSat(this.Sys,this.P, phi, depth, tau, ipts);
        end
        
        
        function [out] = PlotRobustMap(this, phi, params, ranges, delta, options_in)
            % Plot robust satisfaction vs 1 or 2 parameters.
            
            % if zero argument, returns an option structure with defaults
            if nargin==1
                out = struct('contour', 1, 'style',[]);
                return
            end
            
            out = figure;
            
            % no option, use defaults
            if ~exist('options_in','var')
                options_in = struct();
            end
            
            % option provided, make sure all fields are initialized
            options = this.PlotRobustMap();
            opt_in_fields = fieldnames(options_in);
            for  ifld=1:numel(opt_in_fields)
                options.(opt_in_fields{ifld}) = options_in.(opt_in_fields{ifld});
            end
            
            
            switch(nargin)
                case 2
                    this.CheckSpec(phi);
                    figure;
                    SplotProp(this.P, phi, options);
                    return;
                case 4
                    delta = 10;
            end
            
            this.P = CreateParamSet(this.P, params, ranges);
            this.P = Refine(this.P, delta,1);
            this.Sim();
            this.CheckSpec(phi);
            
            SplotProp(this.P, phi, options);
            
        end
        
        %% Experimental
        function report = Analysis(this)
            
            STL_ReadFile('stlib.stl');
            
            % Simple analysis
            % Checks zero signals
            req_zero = STL_Formula('phi', 'req_zero');
            report.res_zero = STL_EvalTemplate(this.Sys, req_zero, this.P, this.P.traj, {'x_'});
            sigs_left = [report.res_zero.some report.res_zero.none];
            
            if ~isempty(report.res_zero.all)
                fprintf('---------------------------------------------\n' )
                fprintf('The following signals were found to be zero: \n' )
                fprintf('---------------------------------------------\n' )
                disp(report.res_zero.all);
            end
            
            % Checks constant signals
            req_stable = STL_Formula('phi', 'req_stable');
            report.res_stable = STL_EvalTemplate(this.Sys, req_stable, this.P, this.P.traj, {'x_'}, sigs_left);
            if ~isempty(report.res_stable.all)
                fprintf('------------------------------------------------\n' )
                fprintf('The following signals were found to be constant: \n' )
                fprintf('------------------------------------------------\n' )
                disp(report.res_stable.all');
            end
            sigs_left = [report.res_stable.some report.res_stable.none];
            
            % Checks non-decreasing signals
            req_inc = STL_Formula('phi', 'req_inc');
            report.res_inc = STL_EvalTemplate(this.Sys, req_inc, this.P, this.P.traj, {'x_'}, sigs_left);
            if ~isempty(report.res_inc.all)
                fprintf('-----------------------------------------------------\n' )
                fprintf('The following signals were found to be non-decreasing: \n' )
                fprintf('-----------------------------------------------------\n' )
                disp(report.res_inc.all');
            end
            
            % Checks non-increasing signals
            req_dec = STL_Formula('phi', 'req_dec');
            report.res_dec = STL_EvalTemplate(this.Sys, req_dec, this.P, this.P.traj, {'x_'}, sigs_left);
            if ~isempty(report.res_dec.all)
                fprintf('-----------------------------------------------------\n' )
                fprintf('The following signals were found to be non-increasing: \n' )
                fprintf('-----------------------------------------------------\n' )
                disp(report.res_dec.all');
            end
            
            % Checks non-increasing signals
            req_dec = STL_Formula('phi', 'req_pwc');
            report.res_dec = STL_EvalTemplate(this.Sys, req_dec, this.P, this.P.traj, {'x_'}, sigs_left);
            if ~isempty(report.res_dec.all)
                fprintf('-------------------------------------------------------\n' )
                fprintf('The following signals were found to be piecewise stable: \n' )
                fprintf('-------------------------------------------------------\n' )
                disp(report.res_dec.all');
            end
            
            % Search for steps
            req_step = STL_Formula('phi', 'req_step');
            report.res_step = STL_EvalTemplate(this.Sys, req_step, this.P, this.P.traj, {'x_step_'}, sigs_left);
            if ~isempty(report.res_step.all)
                fprintf('---------------------------------------------\n' )
                fprintf('The following signals appear to feature steps: \n' )
                fprintf('---------------------------------------------\n' )
                disp(report.res_step.all');
                sigs_steps = report.res_step.all;
            end
            
            % Search for spikes
            req_spike = STL_Formula('phi', 'req_spike');
            report.res_spike = STL_EvalTemplate(this.Sys, req_spike, this.P, this.P.traj, {'x_'}, sigs_left);
            if ~isempty(report.res_spike.all)
                fprintf('---------------------------------------------------------\n' )
                fprintf('The following signals appear to feature spikes or valleys: \n' )
                fprintf('---------------------------------------------------------\n' )
                disp(report.res_spike.all');
                sigs_spikes = report.res_spike.all;
            end
            
            
            % Dual analysis
            % Correlation between steps and spikes
            
            if ~isempty(report.res_spike.all)
                
                req_steps_and_spikes = STL_Formula('phi', 'req_steps_and_spikes');
                report.res_steps_and_spikes = STL_EvalTemplate(this.Sys, req_steps_and_spikes, this.P, this.P.traj, {'x_step_','x_'}, {sigs_steps,sigs_spikes});
                if ~isempty(report.res_steps_and_spikes.all)
                    fprintf('---------------------------------------------------------------------\n' )
                    fprintf('The following pairs of signals appear to be correlated (step => spike): \n' )
                    fprintf('---------------------------------------------------------------------\n' )
                    for i_res= 1:numel(report.res_steps_and_spikes.all)
                        sigs = report.res_steps_and_spikes.all{i_res};
                        disp([ sigs{1} ',' sigs{2}]);
                    end
                end
            end
            
        end
        
        
        
        %% Sensitivity analysis
        function [mu, mustar, sigma] = SensiSpec(this, phi, params, ranges, opt)
            % SensiSpec Sensitivity analysis of a formula to a set of parameters
            this.ResetParamSet();
            opt.tspan = this.Sys.tspan;
            opt.params = FindParam(this.Sys,params);
            opt.lbound = ranges(:,1)';
            opt.ubound = ranges(:,2)';
            opt.plot = 2;
            
            [mu, mustar, sigma, Pr]= SPropSensi(this.Sys, this.P, phi, opt);
            this.P=Pr;
        end
        
        function [monotonicity, Pr, EE] = ChecksMonotony(this, phi, params, ranges, opt)
            % ChecksMonotony performs a quick check to infer monotonicity of a formula wrt parameters
            opt.tspan = this.Sys.tspan;
            opt.params = FindParam(this.Sys,params);
            opt.lbound = ranges(:,1)';
            opt.ubound = ranges(:,2)';
            opt.plot = 0;
            P0 = Sselect(this.P,1);
            
            [~, ~, ~, Pr, EE]= SPropSensi(this.Sys, P0, phi, opt);
            monotonicity = all(EE'>=0)-all(EE'<=0); % 1 if all positive, -1 if all negative, 0 otherwise
        end
        
        %% Printing
        function PrintSpecs(this)
            disp('Specifications:')
            disp('--------------')
            keys = this.Specs.keys;
            for is = 1:numel(keys)
                fprintf('%s\n',keys{is});
            end
            disp(' ');
        end
        
        function PrintAll(this)
            this.PrintParams();
            this.PrintSignals();
            this.PrintSpecs();
        end
        
        function disp(this)
            disp(['BreachSystem with name ' this.Sys.name '.']);
        end
        
        %% GUI
        
        
        function new_phi  = AddSpecGUI(this)
            signals = this.Sys.ParamList(1:this.Sys.DimX);
            new_phi = STL_TemplateGUI('varargin', signals);
            if isa(new_phi, 'STL_Formula')
                this.Specs(get_id(new_phi)) = new_phi;
            end
        end
          
        
        function gui = RunGUI(this)
            P.P = this.P;
            phis=  this.Specs.keys;
            
            Psave(this.Sys, 'Pthis', this.P);
            if (~isempty(phis))
                Propsave(this.Sys, phis{:});
            end
            gui = Breach(this);
            
        end
        
        function ResetFiles(this)
             system(['rm -f ' this.Sys.name '_param_sets.mat']); 
             system(['rm -f ' this.Sys.name '_properties.mat']);
        end
        
    end
end
