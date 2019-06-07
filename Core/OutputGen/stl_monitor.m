classdef stl_monitor < req_monitor
    
    properties
        P0
        P
        Sys
        formula
        predicates
        formula_id
        inout  % input output mode
        relabs % relative absolute mode
        
        rob_map
        diag_map
        formula_names_map 
        
    end
    
    methods
        function this = stl_monitor(formula)
            this.inout = '';
            this.relabs = '';
            if ischar(formula)
                this.formula= STL_Formula(STL_NewID('phi'), formula);
            elseif isa(formula,'STL_Formula')
                this.formula= formula;
            else
                error('stl_monitor:bad_constructor_arg', ...
                    'stl_monitor constructor requires a string or STL_Formula as argument.')
            end
            this.formula_id = get_id(this.formula);
            this.name =  get_id(this.formula);
            
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
        
        function [] = set_mode(this, flag1, flag2)
            switch flag1
                case {'in','out'}
                    this.inout = flag1;
                otherwise
                    this.inout = '';
            end
            switch flag2
                case {'rel','abs'}
                    this.relabs = flag2;
                otherwise
                    this.relabs = '';
            end
        end
        
        function [flag1,flag2] = get_mode(this)
            flag1 = this.inout;
            flag2 = this.relabs;
        end
        
        function [time, Xout] = computeSignals(this, time, X, p)
            this.init_tXp(time,X,p);
            Xout = zeros(numel(this.signals), numel(time));
            
            % compute predicate values
            if ~isempty(this.predicates)
                for ip = 1:numel(this.predicates)
                    [~, Xout(ip,:)] = this.get_standard_rob(this.predicates{ip}, time);
                end
            end
            
            % compute robustnes of top formula
            [time, Xout(end,:)] = this.get_standard_rob(this.formula, time);
            
        end
        
        function [verdict, rob_map, diag_map, formula_names_map] = explain(this,time,X,p)
            
            if nargin>1
                this.init_tXp(time,X,p);
            end
            
            rob_map = containers.Map;
            phi = this.formula;

            [~,~,rob_map] = STL_Eval_IO_Rob(this.Sys, phi, this.P, this.P.traj{1}, 'out', 'rel', rob_map);
            diag_map = containers.Map;
            
            formula_names_map = containers.Map;
            formula_names_map = get_formula_name_map(phi, formula_names_map);
            
            top_signal = rob_map(get_id(phi));
            
            val = top_signal.values(1);
            if(val < 0)
                verdict = 0;
            else
                verdict = 1;
            end
            
            implicant = BreachImplicant;
            implicant = implicant.addInterval(0, 0);
            implicant = implicant.addSignificantSample(0, val);
            
            id = get_id(phi);
            diag_map(id) = implicant;
            
            [this.formula, diag_map] = this.diag(phi, rob_map, diag_map, verdict);
            
            this.rob_map = rob_map;
            this.diag_map = diag_map;
            this.formula_names_map = formula_names_map;
            
        end
        
        function plot_diagnosis(this, F)
            % Assumes F has data about this formula
            [verdict, rob_map, diag_map, formula_names_map] = this.explain();
            
            nb_plots = rob_map.size(1);
            keys = rob_map.keys();

            if (verdict)
                color = 'g';
            else
                color = 'r';
            end
            
            traj = this.P.traj{1}; 
            t0 = traj.time;
            for i=1:nb_plots
                id = keys{i};
                h = F.AddAxes();
                axis(h);
                signal = rob_map(id);
                stairs(signal.times, signal.values);
                
                grid on;               
                formula_name = formula_names_map(id);
                title(formula_name, 'Interpreter', 'none');
                
                ylim = get(h, 'YLim');
                ylim_bot = ylim(1);
                ylim_top = ylim(2);
                
                implicant = diag_map(id);
                size = implicant.getIntervalsSize();
                for j=1:size
                    interval = implicant.getInterval(j);
                    x = interval.begin;
                    y = interval.end;
                    if (x == y)
                        line([x x],[ylim_bot ylim_top],'Color',color);
                    elseif (y > x)
                        p = patch([x y y x], [ylim_bot ylim_bot ylim_top ylim_top], color);
                        alpha(p, 0.05);
                        set(p,'EdgeColor','none');
                    end
                end
                samples = implicant.getSignificantSamples();
                for j=1:length(samples)
                    sample = samples(j);
                    hold on;
                    plot(sample.time, sample.value, 'x');
                end
                
                set(h, 'XLim', [0 t0(end)]);
            end
            
        end
        
        function [v, t, Xout] = eval(this, t, X,p)
            [t, Xout] = this.computeSignals(t, X,p);
            v = Xout(end,1);
        end
        
        function init_tXp(this, t, X, p)
            this.P.traj{1}.time = t;
            
            if isempty(this.signals_in)
                this.P.traj{1}.X = 0*t;
            else
                this.P.traj{1}.X = X;
            end
            
            if nargin>=4&&~isempty(p)
                this.P0 = SetParam(this.P, this.params,p);
            else
                this.P0 = this.P;
            end
        end
        
        function [time, rob] = get_standard_rob(this, phi, time)
            switch this.inout
                case {'in','out'}
                    switch this.relabs
                        case {'rel','abs'}
                            [rob, time] = STL_Eval_IO(this.Sys, phi, this.P0, this.P.traj{1}, this.inout, this.relabs, time);
                        otherwise
                            [rob, time] = STL_Eval_IO(this.Sys, phi, this.P0, this.P.traj{1}, this.inout, 'rel', time);
                    end
                otherwise
                    [rob, time] = STL_Eval(this.Sys, phi, this.P0,this.P.traj{1},time);
            end
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
        function [phi, diag_map] = diag(this, phi, rob, diag_map, flag)
            
            in_implicant = diag_map(get_id(phi));
            samples = in_implicant.getSignificantSamples();
            id = get_id(phi);
            
            psis = get_children(phi);
            switch(get_type(phi))
                case 'predicate'
                    signal_names = STL_ExtractSignals(phi);
                    for i=1:length(signal_names)
                        signal_name = signal_names{i};
                        signal = rob(signal_name);
                        if(~diag_map.isKey(signal_name))
                            out_implicant = BreachImplicant;
                            intervals = in_implicant.getIntervals();
                            for j=1:length(intervals)
                                interval = intervals(j);
                                out_implicant = out_implicant.addInterval(interval.begin, interval.end);
                            end
                            samples = in_implicant.getSignificantSamples();
                            for j=1:length(samples)
                                sample = samples(j);
                                value = interp1(signal.times, signal.values, sample.time, 'previous');
                                out_implicant = out_implicant.addSignificantSample(sample.time, value);
                            end
                            diag_map(signal_name) = out_implicant;
                        end
                    end
                    
                case 'not'
                    signal = rob(get_id(psis{1}));
                    if(flag)
                        [implicant] = BreachDiagnostics.diag_not_t(signal, in_implicant, samples);
                    else
                        [implicant] = BreachDiagnostics.diag_not_f(signal, in_implicant, samples);
                    end
                    diag_map(get_id(psis{1})) = implicant;
                    [psis{1}, diag_map] = this.diag(psis{1}, rob, diag_map, ~flag);
                case 'or'
                    signal1 = rob(get_id(psis{1}));
                    signal2 = rob(get_id(psis{2}));
                    if(flag)
                        [implicant1, implicant2] = BreachDiagnostics.diag_or_t(signal1, signal2, in_implicant, samples);
                    else
                        [implicant1, implicant2] = BreachDiagnostics.diag_or_f(signal1, signal2, in_implicant, samples);
                    end
                    diag_map(get_id(psis{1})) = implicant1;
                    diag_map(get_id(psis{2})) = implicant2;
                    [psis{1}, diag_map] = this.diag(psis{1}, rob, diag_map, flag);
                    [psis{2}, diag_map] = this.diag(psis{2}, rob, diag_map, flag);
                case 'and'
                    signal1 = rob(get_id(psis{1}));
                    signal2 = rob(get_id(psis{2}));
                    if(flag)
                        [implicant1, implicant2] = BreachDiagnostics.diag_and_t(signal1, signal2, in_implicant, samples);
                    else
                        [implicant1, implicant2] = BreachDiagnostics.diag_and_f(signal1, signal2, in_implicant, samples);
                    end
                    diag_map(get_id(psis{1})) = implicant1;
                    diag_map(get_id(psis{2})) = implicant2;
                    [psis{1}, diag_map] = this.diag(psis{1}, rob, diag_map, flag);
                    [psis{2}, diag_map] = this.diag(psis{2}, rob, diag_map, flag);
                case '=>'
                    signal1 = rob(get_id(psis{1}));
                    signal2 = rob(get_id(psis{2}));
                    if(flag)
                        [implicant1, implicant2] = BreachDiagnostics.diag_implies_t(signal1, signal2, in_implicant, samples);
                    else
                        [implicant1, implicant2] = BreachDiagnostics.diag_implies_f(signal1, signal2, in_implicant, samples);
                    end
                    diag_map(get_id(psis{1})) = implicant1;
                    diag_map(get_id(psis{2})) = implicant2;
                    [psis{1}, diag_map] = this.diag(psis{1}, rob, diag_map, flag);
                    [psis{2}, diag_map] = this.diag(psis{2}, rob, diag_map, flag);
                case 'always'
                    signal = rob(get_id(psis{1}));
                    I = this.get_interval();
                    bound.begin = I(1);
                    bound.end = min(I(2),max(signal.times));
                    if(flag)
                        [implicant] = BreachDiagnostics.diag_alw_t(signal, bound, in_implicant, samples);
                    else
                        [implicant] = BreachDiagnostics.diag_alw_f(signal, bound, in_implicant, samples);
                    end
                    
                    diag_map(get_id(psis{1})) = implicant;
                    [psis{1}, diag_map] = this.diag(psis{1}, rob, diag_map, flag);
                case 'eventually'
                    signal = rob(get_id(psis{1}));
                    I = this.get_interval();
                    bound.begin = I(1);
                    bound.end = min(I(2),max(signal.times));
                    if(flag)
                        [implicant] = BreachDiagnostics.diag_ev_t(signal, bound, in_implicant, samples);
                    else
                        [implicant] = BreachDiagnostics.diag_ev_f(signal, bound, in_implicant, samples);
                    end
                    diag_map(get_id(psis{1})) = implicant;
                    [psis{1}, diag_map] = this.diag(psis{1}, rob, diag_map, flag);
            end
        end
        
        function assign_params(this, p)
            if nargin==1
               p = GetParam(this.P, this.params);
            end
            % assign_params fetch parameters and assign them in the current context
            for ip = 1:numel(this.params)
                assignin('caller', this.params{ip},p(ip));
            end
        end
        
        function I = get_interval(this__)
            this__.assign_params();            
            I = eval(get_interval(this__.formula));                
        end
                
        function init_P(this)
            % init_P construct legacy structure from signals and
            % parameters names
            if isempty(this.signals_in)
                sigs = {'dumx'};
            else
                sigs= this.signals_in;
            end
            
            this.Sys = CreateExternSystem([this.formula_id '_Sys'], sigs, this.params, this.p0);
            this.P = CreateParamSet(this.Sys);
            
            traj.param = zeros(1,numel(sigs)+numel(this.params));
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