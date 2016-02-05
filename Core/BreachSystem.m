classdef BreachSystem < BreachSet
    %BreachSystem    A simplified class API for Breach.
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
            end
            [this.P, val] = SEvalProp(this.Sys,this.P,spec);
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
            SplotSat(this.Sys,this.P, phi, depth, tau, ipts )
        end
        
        
        %% Mining
        function [p, rob] = MineSpec(this, phi, falsif_opt, prop_opt, iter_max)
            [p, rob, Pr] = ReqMining(this.Sys, phi, falsif_opt, prop_opt, iter_max);
            this.P = Pr;
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
            disp(['BreachSystem with name ' this.Sys.name '.']);
        end
        
        %% GUI
        function RunGUI(this)
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
