classdef stl_monitor < req_monitor
    % req_monitor wrapper for class STL_Formula
    properties
        P0
        P
        Sys
        formula
        predicates
        formula_id
        inout  % input output mode
        relabs % relative absolute mode
        
        verdict
        rob_map
        diag_map
        formula_names_map 
        
    end
    
    methods
        function this = stl_monitor(formula)
            this.inout = 'out';
            this.relabs = 'rel';
            if ischar(formula)
                this.formula= STL_Formula(STL_NewID('phi'), formula);
            elseif isa(formula,'STL_Formula')
                this.formula= formula;
            elseif isa(formula, 'stl_monitor')
                this.formula = formula.formula;
            else
                error('stl_monitor:bad_constructor_arg', ...
                    'stl_monitor constructor requires a string or STL_Formula as argument.')
            end
            this.formula_id = get_id(this.formula);
            this.name =  get_id(this.formula);
            %this.params_out = {['rho_'  this.name]};
            %this.params_out_domain = BreachDomain();
            
            % collect signals and params names
            [this.signals_in, this.params, this.p0] = STL_ExtractSignals(this.formula);
            input_signals = get_in_signal_names(this.formula);
            output_signals = get_out_signal_names(this.formula);
            if isempty(input_signals)&&isempty(output_signals) % if no IO specified, everybody is output and nobody input
                this.formula = set_out_signal_names(this.formula, this.signals_in);
            end
            this.init_P();            
        end
        
        function set_in_signal_names(this, sigs)
           this.formula = set_in_signal_names(this.formula,sigs,1);
        end
        
        function set_out_signal_names(this, sigs)
           this.formula = set_out_signal_names(this.formula,sigs,1);
        end        
        
        function sigs_in = get_in_signal_names(this)            
            sigs_in = get_in_signal_names(this.formula);
        end
        
        function sigs_out = get_out_signal_names(this)            
           sigs_out = get_out_signal_names(this.formula);
        end
                        
        function status = set_mode(this, flag1, flag2)            
            input_signals = get_in_signal_names(this.formula);
            if ~isempty(input_signals)                                
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
                status = 1;
            else
                status = 0;   % modes are irrelevant
            end
        end
        
        function [flag1,flag2] = get_mode(this)
            flag1 = this.inout;
            flag2 = this.relabs;
        end    
        
        function [time, Xout] = computeSignals(this, time, X, p)
            this.init_tXp(time,X,p);
            
           % compute robustnes of top formula at time 0
           [time, Xout] = this.get_standard_rob(this.formula, 0);
        end
        
        function [v, t, Xout] = eval(this, t, X,p)
            [t, Xout] = this.computeSignals(t, X,p);
            v = Xout(end,1);
        end
        
        function explain(this,time,X,p)
            
            if nargin>1
                this.init_tXp(time,X,p);
            end
            
            this.rob_map = containers.Map;
            phi = this.formula;

            [~,~,this.rob_map] = STL_Eval_IO_Rob(this.Sys, phi, this.P0, this.P.traj{1}, 'out', 'rel', this.rob_map);
            this.diag_map = containers.Map;
            
            this.formula_names_map = containers.Map;
            this.formula_names_map = get_formula_name_map(phi, this.formula_names_map);
            
            top_signal = this.rob_map(get_id(phi));
            
            val = top_signal.values(1);
            if(val < 0)
                this.verdict = 0;
            else
                this.verdict = 1;
            end
            
            implicant = BreachImplicant;
            implicant = implicant.addInterval(0, 0);
            implicant = implicant.addSignificantSample(0, val);
            
            id = get_id(phi);
            this.diag_map(id) = implicant;
            
            [this.formula, this.diag_map] = this.diag(phi, this.rob_map, this.diag_map, this.verdict);
            
        end
        
        function plot_diagnostics(this, F, phi)
            if nargin<3
                phi=this.formula;
            end
            
            this.plot_full_diagnostics(F,phi); 
        end
        
        function plot_full_diagnostics(this,F,phi)
            % Assumes F has data about this formula
            
            itraj = F.BrSet.P.traj_ref(F.ipts);
            time = F.BrSet.P.traj{itraj}.time;
            X = F.BrSet.GetSignalValues(this.signals_in, itraj);
            p = F.BrSet.GetParam(this.params, F.ipts);
            this.explain(time, X,p);
            
            if nargin<3
                phi=this.formula;
            end
            
            subs = STL_Break(phi); 
            plotted = {};
            for is = numel(subs):-1:1
                subphi = subs(is);
                id = get_id(subphi);
                implicant = this.diag_map(id);
                
                if ~isempty(implicant.Intervals)&&~ismember(id, plotted)
                    h = F.AddAxes();
                    this.plot_implicant(h, id);
                    plotted = [plotted id];
                end
            end
            
            for is = 1:numel(this.signals_in)
                id = this.signals_in{is};
                implicant = this.diag_map(id);
                if ~isempty(implicant.Intervals)&&~ismember(id, plotted) 
                    h = F.AddAxes();
                    this.plot_implicant(h, id);                    
                    plotted = [plotted id];
                end
            end
        end
           
        function plot_implicant(this, ax, id)
            
            if (this.verdict)
                color = 'g';
            else
                color = 'r';
            end
            
            axis(ax);
            signal = this.rob_map(id);
            stairs(signal.times, signal.values);
            grid on;
            hold on;
            formula_name = this.formula_names_map(id);
            if ~ismember(id, this.signals_in)
                l = legend([id ': ' formula_name]);
            else
                l = legend(id);
            end
            set(l, 'Interpreter', 'none');
            ylim = get(ax, 'YLim');
            ylim_bot = ylim(1);
            ylim_top = ylim(2);
            
            implicant = this.diag_map(id);
            size = implicant.getIntervalsSize();
            for j=1:size
                interval = implicant.getInterval(j);
                x = interval.begin;
                y = interval.end;
                if (x == y)
                    line([x x],[ylim_bot ylim_top],'Color',color);
                elseif (y > x)
                    p = patch([x y y x], [ylim_bot ylim_bot ylim_top ylim_top], color);
                    alpha(p, 0.3);
                    set(p,'EdgeColor','none');
                end
            end
            samples = implicant.getSignificantSamples();
            %for j=1:length(samples)
            if ~isempty(samples)   
                sample = samples(1);
                hold on;
                plot(sample.time, sample.value, 'x', 'Color',color);            
            end
            t0 = this.P.traj{1}.time;
            set(ax, 'XLim', [t0(1) t0(end)]);
            
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
        
        function [time, rob] = get_standard_rob(this, phi, time, IA_flag)
            
            if nargin<4
                IA_flag = true;
            end
            if IA_flag
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
            else
                [rob, time] = STL_Eval(this.Sys, phi, this.P0,this.P.traj{1},time);
            end
        end
        
        function is_sensitive = ...
                get_structural_sensitivity(this, time, X1, X2, p)
            
            % init for given time
            this.init_tXp(time, X1, p);
            traj1 = this.P.traj{1};
            this.init_tXp(time, X2, p);
            traj2 = this.P.traj{1};
            
            % evaluate sensitivity at time 0
            is_sensitive = ...
                STL_Eval_Structural_Sensitivity(this.Sys, this.formula, ...
                this.P0, traj1, traj2, this.inout, this.relabs, 0);
        end
        
        function varargout = disp(this)
            phi = this.formula;
            st = sprintf(['%s := %s\n'], get_id(phi), disp(phi,1));
            
            if ~strcmp(get_type(phi),'predicate')
                st_pred = [];
                preds = STL_ExtractPredicates(phi);
                for ip = 1:numel(preds)
                    id{ip} = get_id(preds(ip));
                    status(ip)= STL_CheckID(id{ip});
                end
                
                if any(status==1)
                    st_pred = '  where \n';
                    for ip = 1:numel(preds)
                        if status(ip)==1
                            st_pred =   sprintf([ st_pred '%s := %s \n' ],id{ip},disp(preds(ip)));
                        end
                    end
                end
                st = [st st_pred];
            end
     
            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end            
            
        end
    
   
    function assign_params(this, p)
            if nargin==1
               p = GetParam(this.P0, this.params);
            end
            % assign_params fetch parameters and assign them in the current context
            for ip = 1:numel(this.params)
                assignin('caller', this.params{ip},p(ip));
            end
        end
    end
    methods (Access=protected)
        function [phi, diag_map] = diag(this, phi, rob, diag_map, flag)
            
            in_implicant = diag_map(get_id(phi));
            samples = in_implicant.getSignificantSamples();
          
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
                    I = this.get_interval(phi);
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
                    I = this.get_interval(phi);
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
        
        function [mnu, diag_map] = get_mnu(this, phi, rob, diag_map, flag)
            
            in_implicant = diag_map(get_id(phi));
            samples = in_implicant.getSignificantSamples();
          
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
                    I = this.get_interval(phi);
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
                    I = this.get_interval(phi);
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
        
   
        
        function I = get_interval(this__, phi___)
            if nargin==1
                phi___ = this__.formula;
            end
            this__.assign_params();            
            I = eval(get_interval(phi___));                
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