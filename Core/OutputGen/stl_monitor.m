classdef stl_monitor < output_gen
    
    properties
        P
        Sys
        formula
        predicates
        formula_id
    end
    
    methods
        function this = stl_monitor(formula)
            if ischar(formula)
                this.formula= STL_Formula(STL_NewID('phi'), formula);
            elseif isa(formula,'STL_Formula')
                this.formula= formula;
            else
                error('stl_monitor:bad_constructor_arg', ...
                    'stl_monitor constructor requires a string or STL_Formula as argument.')
            end
            this.formula_id = get_id(this.formula);
            
            % collect signals and params names
            [this.signals_in, this.params, this.p0] = STL_ExtractSignals(this.formula);
            
            % Outputs
            if ~strcmp(get_type(this.formula), 'predicate')
                this.signals = {};
                preds = STL_ExtractPredicates(this.formula);
                for ip = 1:numel(preds)
                    if ~STL_CheckID(get_id(preds(ip)))   % predicate does not exist as formula, create it
                        pred = STL_Formula( [get_id(this.formula) '_predicate_' num2str(ip)], preds(ip));
                    else
                        pred = preds(ip);
                    end
                    this.predicates{ip} = pred;
                    this.signals = [this.signals {get_id(this.predicates{ip})}];
                end
                this.signals =  [this.signals {get_id(this.formula)}];
            else
                this.signals = {get_id(this.formula)};
            end
            
            this.init_P();
            
        end
        
        function [tau, Xout] = computeSignals(this, t, X, p, tau)
            if ~exist('tau', 'var')||isempty(tau)
                tau = t;
            end
            this.P.traj{1}.X = X;
            this.P.traj{1}.time = t;
            if nargin>=4&&~isempty(p)
                P0 = SetParam(this.P, this.params,p); % really? P0 gets traj removed,..., gotta get rid of all this non-sense one day
            else
                P0 = this.P;
            end
            Xout = zeros(numel(this.signals), numel(tau));
            Xout(end,:) = STL_Eval(this.Sys, this.formula, P0,this.P.traj{1}, tau);
            if ~isempty(this.predicates)
                for ip = 1:numel(this.predicates)
                    Xout(ip,:) = STL_Eval(this.Sys, this.predicates{ip}, P0,this.P.traj{1}, tau);
                end
            end
        end
        
        function plot_diagnosis(this, F)
            % Assumes F has data about this formula
            if ~isempty(this.predicates)
                F.AddAxes();
                F.AddSignals(this.predicates);
            end
        end
        
        function [v, Xout] = eval(this, t, X,p)
            [~, Xout] = this.computeSignals(t, X,p);
            v = Xout(end,1);
        end
        
        function st = disp(this)
            phi = this.formula;
            st = sprintf(['%s := %s\n'], get_id(phi), disp(phi,1));
            
            if ~strcmp(get_type(phi),'predicate')
                st_pred = [];
                preds = STL_ExtractPredicates(phi);
                for ip = 1:numel(preds)
                    id = get_id(preds(ip));
                    if STL_CheckID(id)
                        if isempty(st_pred)
                            st_pred = '  where \n';
                        end
                        st_pred =   sprintf([ st_pred '%s := %s \n' ],id,disp(preds(ip)));
                    end
                end
                st = [st st_pred];
            end
            
            if nargout == 0
                fprintf(st);
            end
        end
    end
    
    
    methods (Access=protected)
        
        function assign_params(this, p)
            % assign_params fetch parameters and assign them in the current context
            for ip = 1:numel(this.params)
                assignin('caller', this.params{ip},p(ip));
            end
        end
        
        function get_horizon(this__)
            this.assign_params();
            get_interval(this__.formula);
        end
        
        function init_P(this)
            % init_P construct legacy structure from signals and
            % parameters names
            this.Sys = CreateExternSystem([this.formula_id '_Sys'], this.signals_in, this.params, this.p0);
            this.P = CreateParamSet(this.Sys);
            
            traj.param = zeros(1,numel(this.signals_in)+numel(this.params));
            traj.time = [];
            traj.X = [];
            traj.status = 0;
            
            this.P.traj = {traj};
            this.P.traj_ref = 1;
            this.P.traj_to_compute = [];
            
            % Init domains
            for vv =  [this.signals_in this.signals this.params]
                this.domains(vv{1}) = BreachDomain();
            end
        end
    end
end