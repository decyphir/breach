classdef SignalGenerator < BreachSystem
% NOT USED, kept there for reference and eventual adaptation into proper signal_gen classes    
    properties
        signals
        configs
        param_list
        
        type
        time_horizon
        step_size
        value_range
        min_amplitude
        edge_time
        width_range
        time_diff
        period_range
        offset_range
        decay_rate_range
        amplitude
        base_value
        frequency_range
        spike_amplitude_range
        spike_base_range
        spike_width_range
        spike_time
        duty_cycle
        period
        time_range
        min_duration
        omega
        step_amplitude
        step_time
        gain
        ramp_slope
        ramp_start_time
        
        num_edges
        numCp
        baseName
        spikeName
        edgeName
    end
    
    methods
        
        % Constructor - reads a configuration and produces number signals
        % if given two arguments
        function this = SignalGenerator(varargin)
            
            % Initialize parameters
            this.param_list = {'value_range', ...
                'min_amplitude',...
                'edge_time',...
                'width_range',...
                'time_diff',...
                'period_range',...
                'offset_range',...
                'amplitude',...
                'base_value',...
                'frequency_range',...
                'decay_rate_range',...
                'amplitude',...
                'spike_amplitude_range',...
                'spike_base_range',...
                'spike_width_range',...
                'spike_time',...
                'duty_cycle',...
                'period',...
                'time_range',...
                'min_duration',...
                'omega',...
                'step_amplitude',...
                'step_time',...
                'gain',...
                'ramp_slope',...
                'ramp_start_time',...
                };
            
            % Initialize Breach system infrastructure
            sigs = {'x_','x_ref_'};
            params =  this.param_list;
            params =  {params{:}, 'num_edges', 'numCp', 'config'};
            p0 = zeros(1,numel(params));
            p0(:) = NaN;
            simfn = @(Sys, tspan, p) BreachSimWrapper(this, Sys,tspan, p);
            this.Sys =  CreateExternSystem('SignalGen', sigs, params,p0, simfn);
            this.Sys.tspan = 0:.01:20;
            this.P = CreateParamSet(this.Sys);
            this.configs = InitConfigSignals;
            this.Sys.init_fun = @(P) (InitP(this,P));
            
            icfg=1;
            switch nargin
                case 0
                    this = ImportSignalConfig(this, this.configs{icfg});
                case 1
                    icfg = varargin{1};
                    this = ImportSignalConfig(this, this.configs{varargin{1}});
                case 2
                    icfg = varargin{1};
                    this = ImportSignalConfig(this, varargin{1});
                    this = gen_signals(this,varargin{2});
            end
            this.P= SetParam(this.P, 'config',icfg);
        end
        
        %
        function [t, X] = BreachSimWrapper(this, Sys, tspan, p)
            
            icfg = p(end);
            
            try
                cfg = this.configs{icfg};
            catch
                warning('gen_signal:Unknown_Type', 'Unrecognised signal type.');
                t = tspan;
                X = zeros(2,numel(tspan));
                return;
            end
            
            % At this point we have a configuration cfg
            % Because ComputeTraj called InitP before this, p should have
            % values corresponding to cfg, i.e., no NaN, either default or
            % new values
            
            Sys.p =p;
            
            fn = fieldnames(cfg);
            for ifn = 1:numel(fn)
                f = fn{ifn};
                idxf = FindParam(Sys, f);
                if isfield(cfg.(f),'lower')
                    cfg.(f).lower = p(idxf);
                    cfg.(f).upper = p(idxf);
                end
            end
            
            % by default produces signals with 1000 samples
            if numel(tspan) == 2
                tspan = linspace(0, tspan(end), 1000);
            end
            
            cfg.timeHorizon= tspan(end);
            cfg.stepSize = tspan(2)-tspan(1);
            
            this.P.pts(:,:) =0;
            this.ImportSignalConfig(cfg);
            
            
            this.signals={};
            this.gen_signals(1);
            s = this.signals{1};
            
            if size(s,2)==2
                s(1,3) =0;
            end
            
            t = s(:,1)';
            X = s(:,2:end)';
            
        end
        
        function P = InitP(this,P)
            
            nb_pts= size(P.pts,2);           
            for ipts = 1:nb_pts
                p = P.pts(:,ipts);
                icfg = p(end);
                try
                    cfg = this.configs{icfg};
                catch
                    warning('gen_signal:Unknown_Type', 'Unrecognised signal type.');
                    return;
                end
                
                fn = fieldnames(cfg);
                for ifn = 1:numel(fn)
                    f = fn{ifn};
                    if isfield(cfg.(f),'lower')
                        pvalue= GetParam(P,f);
                        ip = FindParam(P,f);
                        if isnan(pvalue)
                            middle_val = (cfg.(f).upper+cfg.(f).lower)/2;
                            P = SetParam(P, f, middle_val );
                        end
                    end
                end
            end
            
        end
        
        % Do it the Breach way - random signals from one signal
        % configuration
        function this = GenRandomSignalsFromCfg(this, cfgs, numbers)
            this.P = [];
            nb_cfgs = numel(cfgs); 
          
            if isscalar(numbers)
                nb = numbers;
                numbers = 0*cfgs+nb;
            end
            
            for icfg = 1:nb_cfgs
               
                cfg = this.configs{cfgs(icfg)};
                fn = fieldnames(cfg);
                params = {};
                ranges = [];
                
                for ifn = 1:numel(fn)
                    f = fn{ifn};
                    if isfield(cfg.(f),'lower')
                        params = {params{:}, f};
                        ranges = [ranges; cfg.(f).lower cfg.(f).upper];
                    end
                end
                
                P = CreateParamSet(this.Sys, params, ranges);
                P = SetParam(P, 'config',cfgs(icfg));
                P = QuasiRefine(P, numbers(icfg));
                this.P = SConcat(this.P, P);
            end
            this.P = ComputeTraj(this.Sys, this.P, this.Sys.tspan);
            
        end
        
        
        
        
        function this = ImportSignalConfig(this,config)
            
            this.type = config.name;
            this.step_size= config.stepSize;
            this.time_horizon = config.timeHorizon;
            
            % parameters with lower and upper bound
            
            for pname = this.param_list
                param_name = pname{1};
                try
                    this.(param_name).lower= config.(param_name).lower;
                    this.(param_name).upper= config.(param_name).upper;
                catch
                end
            end
            
            
            % single valued parameters
            for spname = {'num_edges', ...
                    'numCp',...
                    'baseName',...
                    'spikeName',...
                    'edgeName',...
                    }
                try
                    this.(spname{1})= config.(spname{1});
                catch
                end
                
            end
            
            %  Init default parameters for the default configuration
            fn = fieldnames(config);
            for ifn = 1:numel(fn)
                f = fn{ifn};
                if isfield(config.(f),'lower')
                    middle_val = (config.(f).upper+config.(f).lower)/2;
                    this.P = SetParam(this.P, f, middle_val );
                end
            end
            
            
            
            
        end
        
        function this = gen_signals(this, number)
            
            try
                switch this.type
                    case 'constant'
                        this = gen_constant(this, number);
                    case 'rising_edge'
                        this = gen_edge(this,'rising', number);
                    case 'falling_edge'
                        this = gen_edge(this,'falling',  number);
                    case 'multi_edge_double'
                        this = gen_multi_edge(this,'double', number);
                    case 'multi_edge_boolean'
                        this = gen_multi_edge(this,'boolean', number);
                    case 'monotonic_ascending'
                        this = gen_monotonic(this, 'ascending', number);
                    case 'monotonic_descending'
                        this = gen_monotonic( this,'descending', number);
                    case 'exponentially_increasing'
                        this = gen_exp( this, 'ascending',number);
                    case 'exponentially_decreasing'
                        this = gen_exp( this,'descending', number);
                    case 'sinusoid'
                        this = gen_sinusoid(this,'plain', number);
                    case 'diverging_sinusoid'
                        this = gen_sinusoid(this,'diverging', number);
                    case 'decaying_sinusoid'
                        this = gen_sinusoid(this,'decaying', number);
                    case 'spike'
                        this = gen_spike(this,'up',  number);
                    case 'valley'
                        this = gen_spike( this,'down', number);
                    case 'spike_on_signal'
                        this = gen_spike_on_this(this, number);
                    case 'hunting'
                        this = gen_sinusoid(this,'plain', number);
                    case 'time_restricted_hunting'
                        this = gen_time_restricted_hunting(this, number);
                    case 'hunting_along_signal'
                        this = gen_hunting_along_this(this, number);
                    case 'hunting_after_edge'
                        this = gen_hunting_after_edge(this, number);
                    case 'ramp'
                        this = gen_ramp(this,number);
                    case 'sawtooth'
                        this = gen_sawtooth(this,number);
                    case 'timed_causal'
                        this = gen_timed_causal(this, number);
                    case 'pulse'
                        this = gen_pulse(this,number);
                    case 'pulse_train'
                        this = gen_pulse_train(this,number);
                    case 'overshoot'
                        this = gen_shoot( this,'over', number);
                    case 'undershoot'
                        this = gen_shoot(this,'under',  number);
                    case 'difference'
                        this = gen_difference(this, number);
                    case 'sticking'
                        this = gen_sticking(this, number);
                    case 'underdamped'
                        this = gen_second_order_step_response(this, 'underdamped',number);
                    case 'second_order_step_response'
                        this = gen_second_order_step_response(this,'undamped',number);
                    otherwise
                        warning('gen_signal:Unknown_Type', 'Unrecognised signal type.');
                        this = gen_default(this,number);
                end
            catch
                warning('gen_signal:Config_Error', 'Configuration problem, generating constant signal with value 0.');
                this = gen_default(this,number);
            end
        end
        
        %% Default signal (unrecognized, wrong config etc)
        function this = gen_default(this,number)
            % constant signal has only one parameter, and an upper and lower bound
            times = 0:this.step_size:this.time_horizon;
            for i=1:number
                s = [transpose(times) transpose(0*ones(1,length(times)))];
                this.signals{end+1} = s;
            end
        end
        
        %% Construct constant signal
        function this = gen_constant(this,number)
            % constant signal has only one parameter, and an upper and lower bound
            l = this.value_range.lower;
            u = this.value_range.upper;
            
            vals = l + (u-l)*rand(1,number);
            times = 0:this.step_size:this.time_horizon;
            for i=1:number
                s = [transpose(times) transpose(vals(1, i)*ones(1,length(times)))];
                this.signals{end+1} = s;
            end
        end
        
        %% Construct an edge signal
        function this = gen_edge(this, type,  number)
            % edge signals have three parameters, amplitude, edge_time, and value_range
            % assumption: uv - ua > lv
            la = this.min_amplitude.lower;
            ua = this.min_amplitude.upper;
            lv = this.value_range.lower;
            uv = this.value_range.upper;
            lt = this.edge_time.lower;
            ut = this.edge_time.upper;
            edge_time_this = lt + (ut-lt)*rand(1, number);
            amps = la + (ua-la)*rand(1,number);
            before = lv + (uv - amps - lv).*rand(1,number);
            after = before + amps;
            times = 0:this.step_size:this.time_horizon;
            for i=1:number
                s = before(i)*ones(1,length(times));
                s(times >= edge_time_this(i)) = after(i);
                if (strcmp(type,'falling'))
                    s = (uv - s) + lv;
                end
                this.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        % construct multiple edges this
        function this = gen_multi_edge( this, type, number)
            % multi edge thiss have following parameters:
            % num_edges = how many edges
            % value_range = range of values
            if (strcmp(type, 'double'))
            % making the range of values larger for when upper-lower == 0 
            % (hack) 
                w = 0.1+0.1*max(abs([this.value_range.lower this.value_range.upper]));
                lv = this.value_range.lower-w;
                uv = this.value_range.upper+w;
            end
            nb_edges = this.num_edges;
            times = 0:this.step_size:this.time_horizon;
            for j=1:number
                edgeTimes = this.step_size + (this.time_horizon - 2*this.step_size)*rand(1,nb_edges);
                edgeTimes = sort(edgeTimes);
                edgeTimes = [0 edgeTimes this.time_horizon];
                for i=1:nb_edges+1
                    if (strcmp(type, 'boolean'))
                        if (i==1)
                            s(i) = 0.0;
                        else
                            s(i+1) = 1 - s(i-1);
                        end
                    else
                        s(i) = lv + (uv-lv)*rand(1,1);
                    end
                end
                % the two lines below do piecewise constant interpolation
                p = @(ts,xs,t) xs(min([find(ts-t>0,1)-1 length(xs)]));
                this.signals{end+1} = [transpose(times) transpose(arrayfun(@(t) p(edgeTimes, s, t), times))];
            end
        end
        
        function this = gen_monotonic( this, type, number)
            % monotonic thiss have the following parameters:
            % value_range = what are the values 

            % making the range of values larger for when upper-lower == 0 
            % (hack) 
            w = 0.1+0.1*max(abs([this.value_range.lower this.value_range.upper]));
            l = this.value_range.lower-w;
            u = this.value_range.upper+w;
            times = 0:this.step_size:this.time_horizon;
            for i=1:number
                s = l + (u-l)*rand(1, length(times));
                if (strcmp(type, 'ascending'))
                    s = sort(s);
                else
                    s = sort(s, 'descend');
                end
                this.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        function this = gen_exp( this, type, number)
            % exponential thiss have following parameters:
            % decay_rate_range = range for decay rates
            % value_range = range for values
            % amplitude = min diff between x(0) and x(inf)
            
            lv = this.value_range.lower;
            uv = this.value_range.upper;
            la = this.min_amplitude.lower;
            ua = this.min_amplitude.upper;
            lr = this.decay_rate_range.lower;
            ur = this.decay_rate_range.upper;
            
            if (min(lr,ur) < 0)
                warning('Decay rate should be positive! Taking absolute values.');
                lr = min(abs(lr),abs(ur));
                ur = max(abs(lr),abs(ur));
            end
            times = 0:this.step_size:this.time_horizon;
            amps = la + (ua-la)*rand(1,number);
            rates = lr + (ur-lr)*rand(1,number);
            before = lv + (uv - amps -lv);
            after = before + amps;
            for i=1:number
                s = after(i) - (after(i)-before(i))*exp(-rates(i).*times);
                if (strcmp(type,'descending'))
                    s = lv + (uv - s);
                end
                this.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        function this = gen_sinusoid(this, type, number)
            % sinusoids have the following parameters:
            % decay_rate_range
            % frequency_range
            lv = this.amplitude.lower;
            uv = this.amplitude.upper;
            lb = this.base_value.lower;
            ub = this.base_value.upper;
            
            rates = zeros(1,number);
            amps = lv + (uv-lv)*rand(1,number);
            bases = lb + (ub-lb)*rand(1,number);
            
            if (~strcmp(type,'plain'))
                lr = this.decay_rate_range.lower;
                ur = this.decay_rate_range.upper;
                if (min(lr,ur) < 0)
                    warning('Decay rate should be positive! Taking absolute values.');
                    lr = min(abs(lr),abs(ur));
                    ur = max(abs(lr),abs(ur));
                end
                rates = lr + (ur - lr)*rand(1,number);
                if (strcmp(type,'decaying'))
                    mult = -1.0;
                else
                    mult = 1.0;
                end
            else
                mult = 0.0;
            end
            lf = this.frequency_range.lower;
            uf = this.frequency_range.upper;
            freqs = lf + (uf-lf)*rand(1,number);
            times = 0:this.step_size:this.time_horizon;
            
            for i=1:number
                s = bases(i) + amps(i)*exp(mult*rates(i)*times).*sin(freqs(i)*times);
                this.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        function this = gen_spike (this,type,number)
            % spike thiss have following parameters:
            % spike_amplitude_range
            % spike_base_range
            % spike_width_range: range of spike_width at 37% of the spike_amplitude
            % spike_time : range of time when spike can happen
            la = this.spike_amplitude_range.lower;
            ua = this.spike_amplitude_range.upper;
            lb = this.spike_base_range.lower;
            ub = this.spike_base_range.upper;
            lw = this.spike_width_range.lower;
            uw = this.spike_width_range.upper;
            lt = this.spike_time.lower;
            ut = this.spike_time.upper;
            
            as = la + (ua-la)*rand(1,number);
            bs = lb + (ub-lb)*rand(1,number);
            ws = lw + (uw-lw)*rand(1,number);
            ts = lt + (ut-lt)*rand(1,number);
            
            times = 0:this.step_size:this.time_horizon;
            if (strcmp(type, 'up'))
                mult = 1;
            else
                mult = -1;
            end
            for i=1:number
                s = bs(i) + mult*as(i)*exp(-((ts(i) - times).^2) / (ws(i)/2));
                this.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        function this = gen_spike_on_this(this, number)
            basethis = copy(this);
            basethis.type = this.baseName;
            basethis = basethis.gen_signals(number);
            spikethis = copy(this);
            spikethis.type = this.spikeName;
            spikethis = spikethis.gen_signals(number);
            for i=1:number
                baseSig = basethis.signals{i};
                spikeSig = spikethis.signals{i};
                this.signals{end+1} = [baseSig(:,1) baseSig(:,2)+spikeSig(:,2)];
            end
        end
        
        function this = gen_time_restricted_hunting(this, number)
            huntingthis = copy(this);
            huntingthis.type = 'hunting';
            huntingthis = gen_signals(huntingthis, number);
            lt = this.time_range.lower;
            ut = this.time_range.upper;
            ld = this.min_duration.lower;
            ud = this.min_duration.upper;
            
            durations = ld + (ud-ld)*rand(1,number);
            lrange = lt + (ut - durations).*rand(1,number);
            urange = lrange + durations;
            for i=1:number
                sig = huntingthis.signals{i};
                sig(sig(:,1)<=lrange(i) | sig(:,1) >=urange(i), 2) = mean(sig(:,2));
                this.signals{end+1} = sig;
            end
        end
        
        function this = gen_hunting_along_this(this, number)
            basethis = copy(this);
            basethis.type = this.baseName;
            basethis = basethis.gen_signals(number);
            huntingthis = gen_time_restricted_hunting(this, number);
            for i=1:number
                baseSig = basethis.signals{i};
                huntingSig = huntingthis.signals{i};
                this.signals{end+1} = [baseSig(:,1) baseSig(:,2)+huntingSig(:,2)];
            end
        end
        
        function this = gen_hunting_after_edge(this, number)
            edgethis = copy(this);
            edgethis.type = this.edgeName;
            edgethis = gen_signals(edgethis, number);
            huntingthis = gen_sinusoid(this,'plain', number);
            for i=1:number
                edgeSig = edgethis.signals{i};
                huntingSig = huntingthis.signals{i};
                z = find(abs(diff(edgeSig(:,2)))>0);
                edgeSig(z:end,2) = edgeSig(z:end,2) + huntingSig(z:end,2);
                this.signals{end+1} = edgeSig;
            end
        end
        
        function this = gen_shoot( this,type, number)
            if (strcmp(type,'over'))
                edgethis = copy(this);
                edgethis = gen_edge(edgethis, 'rising', number);
                mult = 1;
            else
                edgethis = copy(this);
                edgethis = gen_edge(edgethis, 'falling',number);
                mult = -1;
            end
            
            for i=1:number
                u = edgethis.signals{i};
                uSig = u(:,2);
                uTimes = u(:,1);
                edgeTimeIndex = find(abs(diff(uSig))>0, 1, 'first');
                edgeTime = uTimes(edgeTimeIndex);
                stepValue = mult*uSig(edgeTimeIndex+1,1);
                
                if isfield(this,'gain')
                    gain_upper = this.gain.upper*stepValue;
                    gain_lower = this.gain.lower*stepValue;
                else
                    gain_upper = stepValue;
                    gain_lower = stepValue;
                end
                
                if isfield(this,'zeta')
                    transient_cfg = struct(...
                        'stepSize', this.step_size, ...
                        'timeHorizon', this.time_horizon, ...
                        'name', 'underdamped', ...
                        'gain', struct('lower', gain_lower, 'upper', gain_upper), ...
                        'step_time', struct('lower', edgeTime, 'upper', edgeTime),...
                        'step_amplitude', struct('lower', 1, 'upper', 1), ...
                        'omega', struct('lower', this.omega.lower, ...
                        'upper', this.omega.upper),...
                        'zeta', struct('lower', this.zeta.lower, ...
                        'upper', this.zeta.upper));
                else
                    transient_cfg = struct(...
                        'stepSize', this.step_size, ...
                        'timeHorizon', this.time_horizon, ...
                        'name', 'underdamped', ...
                        'gain', struct('lower', gain_lower, 'upper', gain_upper), ...
                        'step_time', struct('lower', edgeTime, 'upper', edgeTime),...
                        'step_amplitude', struct('lower', 1, 'upper', 1), ...
                        'omega', struct('lower', this.omega.lower, ...
                        'upper', this.omega.upper));
                end
                
                transientthis = SignalGenerator();
                transientthis.ImportSignalConfig(transient_cfg);
                transientthis.gen_signals(1);
                
                y = transientthis.signals{1};
                this.signals{end+1} = [uTimes y(:,2) uSig];
            end
        end
        
        function this = gen_difference(this, number)
            % generate two random thiss with sup-norm between them at most a
            % user-specified amount
            ldiff = this.sup_diff.lower;
            udiff = this.sup_diff.upper;
            numCp = this.num_control_points;
            interpType = this.interpolationType;
            lv = this.value_range.lower;
            uv = this.value_range.upper;
            for i=1:number
                cpValues = lv + (uv-lv)*rand(1,numCp);
                cpTimes = linspace(0,this.time_horizon,numCp);
                sigTimes = 0:this.step_size:this.time_horizon;
                if (strcmp(interpType, 'pconst'))
                    p = @(ts,xs,t) xs(min([find(ts-t>0,1)-1 length(xs)]));
                    sig1 = arrayfun(@(t) p(cpTimes, s, t), sigTimes);
                else
                    sig1 = interp1(cpTimes, cpValues, sigTimes, interpType);
                end
                incDec = -1 + 2*randi([0 1], 1, length(sigTimes));
                diffs = incDec.*(ldiff + (udiff-ldiff)*rand(1,length(sigTimes)));
                sig2 = sig1 + diffs;
                this.signals{end+1} = [tranpose(sigTimes) tranpose(sig1) transpose(sig2)];
            end
        end
        
        function this = gen_sticking(this, number)
            constthis = gen_constant(this,number);
            l = this.value_range.lower;
            u = this.value_range.upper;
            lr = this.decay_rate_range.lower;
            ur = this.decay_rate_range.upper;
            rates = lr + (ur-lr)*rand(1,number);
            vals = l + (u-l)*rand(1,number);
            
            for i=1:number
                sig1 = constthis.signals{end+1};
                constSig = sig1(:,2);
                constVal = constSig(1,1);
                sigTimes = sig1(:,1);
                val = vals(i);
                if (val < constVal)
                    % ascending
                    sig2 = constVal - (constVal - val)*exp(-rates(i)*sigTimes);
                else
                    % descending
                    sig2 = (val-constVal)*exp(-rates(i)*sigTimes) + constVal;
                end
                this.signals{end+1} = [sigTimes constSig sig2];
            end
        end
        
        function this = gen_ramp(this, number)
            % assumes slope != 0
            l = this.value_range.lower;
            u = this.value_range.upper;
            ls = this.ramp_slope.lower;
            us = this.ramp_slope.upper;
            lt = this.ramp_start_time.lower;
            ut = this.ramp_start_time.upper;
            la = this.min_amplitude.lower;
            ua = this.min_amplitude.upper;
            slopes = ls + (us-ls)*rand(1, number);
            startTimes = lt + (ut-lt)*rand(1,number);
            amps = la + (ua-la)*rand(1,number);
            
            lvals = l + (u-amps-l).*rand(1,number);
            uvals = lvals + amps;
            times = 0:this.step_size:this.time_horizon;
            % given m, xp, yp, computes y for a given x
            lineY = @(m,xp,yp,x) (m*(x-xp) + yp);
            % given m, xp, yp, computes x for a given y
            lineX = @(m,xp,yp,y) iif(m==0, 0, (((y-yp)/m) + xp));
            for i=1:number
                startValue = lvals(i);
                finalValue = uvals(i);
                if (slopes(i) < 0)
                    startValue = uvals(i);
                    finalValue = lvals(i);
                end
                startTime = startTimes(i);
                finalTime = lineX(slopes(i), startTime, startValue, finalValue);
                givenRamp = @(t) lineY(slopes(i), startTime, startValue, t);
                sig = givenRamp(times);
                sig(times <= startTime) = startValue;
                if (times(end) > finalTime)
                    sig(times >= finalTime) = finalValue;
                end
                this.signals{end+1} = [transpose(times) transpose(sig)];
            end
        end
        
        function this = gen_sawtooth(this,number)
            l = this.value_range.lower;
            u = this.value_range.upper;
            lam = this.min_amplitude.lower;
            uam = this.min_amplitude.upper;
            lw = this.width_range.lower;
            uw = this.width_range.upper;
            lp = this.period_range.lower;
            up = this.period_range.upper;
            la = this.offset_range.lower;
            ua = this.offset_range.upper;
            amps = lam + (uam-lam)*rand(1,number);
            lvals = l + (u - uam - l).*rand(1,number);
            uvals = lvals + amps;
            times = 0:this.step_size:this.time_horizon;
            offsets = la + (ua-la)*rand(1,number);
            widths = lw + (uw-lw)*rand(1,number); % 0<w<1
            periods = lp + (up-lp)*rand(1,number);
            sawt = @(l,u,o,p,w,t) (l + (u-l)*sawtooth(2*pi*(t - o)/p,w));
            for i=1:number
                sig = sawt(lvals(i), uvals(i), offsets(i), ...
                    periods(i), widths(i), times);
                this.signals{end+1} = [transpose(times) transpose(sig)];
            end
        end
        
        function this = gen_timed_causal(this, number)
            % alw( flag1  => ev_[0,t] flag2)
            lt = this.time_diff.lower;
            ut = this.time_diff.upper;
            cpArray = linspace(0, this.time_horizon, this.numCp);
            times = 0:this.step_size:this.time_horizon;
            p = @(ts,xs,t) xs(min([find(ts-t>0,1)-1 length(xs)]));
            for i=1:number
                if (rand(1,1) < 0.01)
                    % 1 % chance of producing a this in which antecedent is never
                    % true
                    sig1cp = zeros(1,length(cpArray));
                    sig2cp = randi([0 1],1,length(cpArray));
                    sig2 = arrayfun(@(t) p(cpArray, sig2cp, t), times);
                    sig1 = arrayfun(@(t) p(cpArray, sig1cp, t), times);
                else
                    sig1cp = randi([0 1], 1, length(cpArray));
                    sig1 = arrayfun(@(t) p(cpArray, sig1cp, t), times);
                    sig2 = zeros(1,length(times));
                    endTime = this.time_horizon - ut;
                    lastInd = find(times < endTime, 1, 'last');
                    matches = find(sig1(1,1:lastInd));
                    dt = lt + (ut-lt)*rand(1,length(matches));
                    for j=1:length(matches)
                        sig1Index = matches(j);
                        sig2Index = find(times < (times(sig1Index) + dt(j)), 1, 'last');
                        % stay high till next cp. so, first find next cp
                        sig2Last = ceil (ceil((sig2Index / length(times))*this.numCp) * ...
                            (length(times)/this.numCp));
                        sig2(sig2Index:sig2Last) = 1;
                    end
                    sig1(lastInd+1:end) = 0;
                end
                this.signals{end+1} = [transpose(times) transpose(sig2(1:numel(times))) transpose(sig1(1:numel(times))) ];
            end
        end
        
        function this = gen_pulse(this, number)
            % assume: ut + ud < time_horizon
            lb = this.base_value.lower;
            ub = this.base_value.upper;
            
            la = this.amplitude.lower;
            ua = this.amplitude.upper;
            
            ld = this.duration.lower;
            ud = this.duration.upper;
            
            lt = this.start_time.lower;
            ut = this.start_time.upper;
            
            times = 0:this.step_size:this.time_horizon;
            startTimes = lt + (ut-lt)*rand(1,number);
            durations = ld + (ud-ld)*rand(1,number);
            bases = lb + (ub-lb)*rand(1,number);
            amplitudes = la + (ua-la)*rand(1,number);
            
            for i=1:number
                sig = bases(i)*ones(1,length(times));
                sig( times > startTimes(i) & ...
                    times < startTimes(i) + durations(i) ) = ...
                    bases(i) + amplitudes(i);
                this.signals{end+1} = [transpose(times) transpose(sig)];
            end
        end
        
        function this = gen_pulse_train(this, number)
            % assume: ut + ud < time_horizon
            lb = this.base_value.lower;
            ub = this.base_value.upper;
            
            la = this.amplitude.lower;
            ua = this.amplitude.upper;
            
            ld = this.duty_cycle.lower;
            ud = this.duty_cycle.upper;
            
            lp = this.period.lower;
            up = this.period.upper;
            
            times = 0:this.step_size:this.time_horizon;
            periods = lp + (up-lp)*rand(1,number);
            dutyCycles = ld + (ud-ld)*rand(1,number);
            bases = lb + (ub-lb)*rand(1,number);
            amplitudes = la + (ua-la)*rand(1,number);
            
            unshifted_pulse = @(T,d,t) sin(((2*pi)/T)*t) > cos(pi*d);
            advance = @(T,d) T*(0.25 - 0.5*d);
            pulse = @(T,d,t) unshifted_pulse(T,d, t + advance(T,d));
            
            for i=1:number
                sig = bases(i) + amplitudes(i)*pulse(periods(i), dutyCycles(i), times);
                this.signals{end+1} = [transpose(times) transpose(sig)];
            end
        end
        
        function this = gen_second_order_step_response(this,type, number)
            la = this.step_amplitude.lower;
            ua = this.step_amplitude.upper;
            ls = this.step_time.lower;
            us = this.step_time.upper;
            lw = this.omega.lower;
            uw = this.omega.upper;
            lg = this.gain.lower;
            ug = this.gain.upper;
            omegas = lw + (uw-lw)*rand(1,number);
            gains = lg + (ug-lg)*rand(1,number);
            amps = la + (ua-la)*rand(1,number);
            sts = ls + (us-ls)*rand(1,number);
            if isfield(this,'zeta')
                zetal = this.zeta.lower;
                zetau = this.zeta.upper;
                zetas = zetal + (zetau-zetal)*rand(1,number);
            else
                switch type
                    case 'underdamped'
                        zetas = 0.2+0.6*rand(1,number);
                    case 'overdamped'
                        zetas = 1.5 + 0.5*rand(1,number);
                    case 'critically_damped'
                        zetas = ones(1,number);
                    case 'undamped'
                        zetas = zeros(1,number);
                end
            end
            
            times = 0:this.step_size:this.time_horizon;
            for i=1:number
                sys = tf([gains(i)*omegas(i)^2], [1 2*zetas(i)*omegas(i) omegas(i)^2]);
                u = amps(i)*heaviside(times-sts(i));
                this.signals{end+1} = [transpose(times) lsim(sys,u,times)];
            end
        end
        
        function plot_signals(this)
            plot_signals(this.signals, this.type);
        end
        
        % Make a copy of a handle object - works because no property is
        % itself a handle object...
        function new = copy(this)
            % Instantiate new object of the same class.
            new = feval(class(this));
            
            % Copy all non-hidden properties.
            p = fieldnames(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end
        
        function disp(this)
            
            disp('Signal Generator with Configurations:')
            disp('-------------------------------------')
            for ip = 1:numel(this.configs)
                fprintf('%i: %s\n', ip, this.configs{ip}.name);
            end
            
        end
        
    end
    
end


function res = iif(b, t, e)
if (b)
    res = t;
else
    res = e;
end
end

