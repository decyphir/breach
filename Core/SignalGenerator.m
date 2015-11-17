classdef SignalGenerator < BreachObject
    
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
        function SG = SignalGenerator(varargin)
            
            % Initialize parameters
            SG.param_list = {'value_range', ...
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
            sigs = {'x','xref'};
            params =  SG.param_list;
            params =  {params{:}, 'num_edges', 'numCp', 'config'};
            p0 = zeros(1,numel(params));
            p0(:) = NaN;
            simfn = @(Sys, tspan, p) BreachSimWrapper(SG, Sys,tspan, p);
            SG.Sys =  CreateExternSystem('SignalGen', sigs, params,p0, simfn);
            SG.Sys.tspan = 0:.01:20;
            SG.P = CreateParamSet(SG.Sys);
            SG.configs = InitConfigSignals;
            SG.Sys.init_fun = @(P) (InitP(SG,P));
            
            icfg=1;
            switch nargin
                case 0
                    SG = ImportSignalConfig(SG, SG.configs{icfg});
                case 1
                    icfg = varargin{1};
                    SG = ImportSignalConfig(SG, SG.configs{varargin{1}});
                case 2
                    icfg = varargin{1};
                    SG = ImportSignalConfig(SG, varargin{1});
                    SG = gen_signals(SG,varargin{2});
            end
            SG.P= SetParam(SG.P, 'config',icfg);
        end
        
        %
        function [t, X] = BreachSimWrapper(SG, Sys, tspan, p)
            
            icfg = p(end);
            
            try
                cfg = SG.configs{icfg};
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
            
            SG.P.pts(:,:) =0;
            SG.ImportSignalConfig(cfg);
            
            
            SG.signals={};
            SG.gen_signals(1);
            s = SG.signals{1};
            
            if size(s,2)==2
                s(1,3) =0;
            end
            
            t = s(:,1)';
            X = s(:,2:end)';
            
        end
        
        function P = InitP(SG,P)
            
            nb_pts= size(P.pts,2);           
            for ipts = 1:nb_pts
                p = P.pts(:,ipts);
                icfg = p(end);
                try
                    cfg = SG.configs{icfg};
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
        function SG = GenRandomSignalsFromCfg(SG, cfgs, numbers)
            SG.P = [];
            nb_cfgs = numel(cfgs); 
          
            if isscalar(numbers)
                nb = numbers;
                numbers = 0*cfgs+nb;
            end
            
            for icfg = 1:nb_cfgs
               
                cfg = SG.configs{cfgs(icfg)};
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
                
                P = CreateParamSet(SG.Sys, params, ranges);
                P = SetParam(P, 'config',cfgs(icfg));
                P = QuasiRefine(P, numbers(icfg));
                SG.P = SConcat(SG.P, P);
            end
            SG.P = ComputeTraj(SG.Sys, SG.P, SG.Sys.tspan);
            
        end
        
        
        
        
        function SG = ImportSignalConfig(SG,config)
            
            SG.type = config.name;
            SG.step_size= config.stepSize;
            SG.time_horizon = config.timeHorizon;
            
            % parameters with lower and upper bound
            
            for pname = SG.param_list
                param_name = pname{1};
                try
                    SG.(param_name).lower= config.(param_name).lower;
                    SG.(param_name).upper= config.(param_name).upper;
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
                    SG.(spname{1})= config.(spname{1});
                catch
                end
                
            end
            
            %  Init default parameters for the default configuration
            fn = fieldnames(config);
            for ifn = 1:numel(fn)
                f = fn{ifn};
                if isfield(config.(f),'lower')
                    middle_val = (config.(f).upper+config.(f).lower)/2;
                    SG.P = SetParam(SG.P, f, middle_val );
                end
            end
            
            
            
            
        end
        
        function SG = gen_signals(SG, number)
            
            try
                switch SG.type
                    case 'constant'
                        SG = gen_constant(SG, number);
                    case 'rising_edge'
                        SG = gen_edge(SG,'rising', number);
                    case 'falling_edge'
                        SG = gen_edge(SG,'falling',  number);
                    case 'multi_edge_double'
                        SG = gen_multi_edge(SG,'double', number);
                    case 'multi_edge_boolean'
                        SG = gen_multi_edge(SG,'boolean', number);
                    case 'monotonic_ascending'
                        SG = gen_monotonic(SG, 'ascending', number);
                    case 'monotonic_descending'
                        SG = gen_monotonic( SG,'descending', number);
                    case 'exponentially_increasing'
                        SG = gen_exp( SG, 'ascending',number);
                    case 'exponentially_decreasing'
                        SG = gen_exp( SG,'descending', number);
                    case 'sinusoid'
                        SG = gen_sinusoid(SG,'plain', number);
                    case 'diverging_sinusoid'
                        SG = gen_sinusoid(SG,'diverging', number);
                    case 'decaying_sinusoid'
                        SG = gen_sinusoid(SG,'decaying', number);
                    case 'spike'
                        SG = gen_spike(SG,'up',  number);
                    case 'valley'
                        SG = gen_spike( SG,'down', number);
                    case 'spike_on_signal'
                        SG = gen_spike_on_SG(SG, number);
                    case 'hunting'
                        SG = gen_sinusoid(SG,'plain', number);
                    case 'time_restricted_hunting'
                        SG = gen_time_restricted_hunting(SG, number);
                    case 'hunting_along_signal'
                        SG = gen_hunting_along_SG(SG, number);
                    case 'hunting_after_edge'
                        SG = gen_hunting_after_edge(SG, number);
                    case 'ramp'
                        SG = gen_ramp(SG,number);
                    case 'sawtooth'
                        SG = gen_sawtooth(SG,number);
                    case 'timed_causal'
                        SG = gen_timed_causal(SG, number);
                    case 'pulse'
                        SG = gen_pulse(SG,number);
                    case 'pulse_train'
                        SG = gen_pulse_train(SG,number);
                    case 'overshoot'
                        SG = gen_shoot( SG,'over', number);
                    case 'undershoot'
                        SG = gen_shoot(SG,'under',  number);
                    case 'difference'
                        SG = gen_difference(SG, number);
                    case 'sticking'
                        SG = gen_sticking(SG, number);
                    case 'underdamped'
                        SG = gen_second_order_step_response(SG, 'underdamped',number);
                    case 'second_order_step_response'
                        SG = gen_second_order_step_response(SG,'undamped',number);
                    otherwise
                        warning('gen_signal:Unknown_Type', 'Unrecognised signal type.');
                        SG = gen_default(SG,number);
                end
            catch
                warning('gen_signal:Config_Error', 'Configuration problem, generating constant signal with value 0.');
                SG = gen_default(SG,number);
            end
        end
        
        %% Default signal (unrecognized, wrong config etc)
        function SG = gen_default(SG,number)
            % constant signal has only one parameter, and an upper and lower bound
            times = 0:SG.step_size:SG.time_horizon;
            for i=1:number
                s = [transpose(times) transpose(0*ones(1,length(times)))];
                SG.signals{end+1} = s;
            end
        end
        
        %% Construct constant signal
        function SG = gen_constant(SG,number)
            % constant signal has only one parameter, and an upper and lower bound
            l = SG.value_range.lower;
            u = SG.value_range.upper;
            
            vals = l + (u-l)*rand(1,number);
            times = 0:SG.step_size:SG.time_horizon;
            for i=1:number
                s = [transpose(times) transpose(vals(1, i)*ones(1,length(times)))];
                SG.signals{end+1} = s;
            end
        end
        
        %% Construct an edge signal
        function SG = gen_edge(SG, type,  number)
            % edge signals have three parameters, amplitude, edge_time, and value_range
            % assumption: uv - ua > lv
            la = SG.min_amplitude.lower;
            ua = SG.min_amplitude.upper;
            lv = SG.value_range.lower;
            uv = SG.value_range.upper;
            lt = SG.edge_time.lower;
            ut = SG.edge_time.upper;
            edge_time_this = lt + (ut-lt)*rand(1, number);
            amps = la + (ua-la)*rand(1,number);
            before = lv + (uv - amps - lv).*rand(1,number);
            after = before + amps;
            times = 0:SG.step_size:SG.time_horizon;
            for i=1:number
                s = before(i)*ones(1,length(times));
                s(times >= edge_time_this(i)) = after(i);
                if (strcmp(type,'falling'))
                    s = (uv - s) + lv;
                end
                SG.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        % construct multiple edges SG
        function SG = gen_multi_edge( SG, type, number)
            % multi edge SGs have following parameters:
            % num_edges = how many edges
            % value_range = range of values
            if (strcmp(type, 'double'))
            % making the range of values larger for when upper-lower == 0 
            % (hack) 
                w = 0.1+0.1*max(abs([SG.value_range.lower SG.value_range.upper]));
                lv = SG.value_range.lower-w;
                uv = SG.value_range.upper+w;
            end
            nb_edges = SG.num_edges;
            times = 0:SG.step_size:SG.time_horizon;
            for j=1:number
                edgeTimes = SG.step_size + (SG.time_horizon - 2*SG.step_size)*rand(1,nb_edges);
                edgeTimes = sort(edgeTimes);
                edgeTimes = [0 edgeTimes SG.time_horizon];
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
                SG.signals{end+1} = [transpose(times) transpose(arrayfun(@(t) p(edgeTimes, s, t), times))];
            end
        end
        
        function SG = gen_monotonic( SG, type, number)
            % monotonic SGs have the following parameters:
            % value_range = what are the values 

            % making the range of values larger for when upper-lower == 0 
            % (hack) 
            w = 0.1+0.1*max(abs([SG.value_range.lower SG.value_range.upper]));
            l = SG.value_range.lower-w;
            u = SG.value_range.upper+w;
            times = 0:SG.step_size:SG.time_horizon;
            for i=1:number
                s = l + (u-l)*rand(1, length(times));
                if (strcmp(type, 'ascending'))
                    s = sort(s);
                else
                    s = sort(s, 'descend');
                end
                SG.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        function SG = gen_exp( SG, type, number)
            % exponential SGs have following parameters:
            % decay_rate_range = range for decay rates
            % value_range = range for values
            % amplitude = min diff between x(0) and x(inf)
            
            lv = SG.value_range.lower;
            uv = SG.value_range.upper;
            la = SG.min_amplitude.lower;
            ua = SG.min_amplitude.upper;
            lr = SG.decay_rate_range.lower;
            ur = SG.decay_rate_range.upper;
            
            if (min(lr,ur) < 0)
                warning('Decay rate should be positive! Taking absolute values.');
                lr = min(abs(lr),abs(ur));
                ur = max(abs(lr),abs(ur));
            end
            times = 0:SG.step_size:SG.time_horizon;
            amps = la + (ua-la)*rand(1,number);
            rates = lr + (ur-lr)*rand(1,number);
            before = lv + (uv - amps -lv);
            after = before + amps;
            for i=1:number
                s = after(i) - (after(i)-before(i))*exp(-rates(i).*times);
                if (strcmp(type,'descending'))
                    s = lv + (uv - s);
                end
                SG.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        function SG = gen_sinusoid(SG, type, number)
            % sinusoids have the following parameters:
            % decay_rate_range
            % frequency_range
            lv = SG.amplitude.lower;
            uv = SG.amplitude.upper;
            lb = SG.base_value.lower;
            ub = SG.base_value.upper;
            
            rates = zeros(1,number);
            amps = lv + (uv-lv)*rand(1,number);
            bases = lb + (ub-lb)*rand(1,number);
            
            if (~strcmp(type,'plain'))
                lr = SG.decay_rate_range.lower;
                ur = SG.decay_rate_range.upper;
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
            lf = SG.frequency_range.lower;
            uf = SG.frequency_range.upper;
            freqs = lf + (uf-lf)*rand(1,number);
            times = 0:SG.step_size:SG.time_horizon;
            
            for i=1:number
                s = bases(i) + amps(i)*exp(mult*rates(i)*times).*sin(freqs(i)*times);
                SG.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        function SG = gen_spike (SG,type,number)
            % spike SGs have following parameters:
            % spike_amplitude_range
            % spike_base_range
            % spike_width_range: range of spike_width at 37% of the spike_amplitude
            % spike_time : range of time when spike can happen
            la = SG.spike_amplitude_range.lower;
            ua = SG.spike_amplitude_range.upper;
            lb = SG.spike_base_range.lower;
            ub = SG.spike_base_range.upper;
            lw = SG.spike_width_range.lower;
            uw = SG.spike_width_range.upper;
            lt = SG.spike_time.lower;
            ut = SG.spike_time.upper;
            
            as = la + (ua-la)*rand(1,number);
            bs = lb + (ub-lb)*rand(1,number);
            ws = lw + (uw-lw)*rand(1,number);
            ts = lt + (ut-lt)*rand(1,number);
            
            times = 0:SG.step_size:SG.time_horizon;
            if (strcmp(type, 'up'))
                mult = 1;
            else
                mult = -1;
            end
            for i=1:number
                s = bs(i) + mult*as(i)*exp(-((ts(i) - times).^2) / (ws(i)/2));
                SG.signals{end+1} = [transpose(times) transpose(s)];
            end
        end
        
        function SG = gen_spike_on_SG(SG, number)
            baseSG = copy(SG);
            baseSG.type = SG.baseName;
            baseSG = baseSG.gen_signals(number);
            spikeSG = copy(SG);
            spikeSG.type = SG.spikeName;
            spikeSG = spikeSG.gen_signals(number);
            for i=1:number
                baseSig = baseSG.signals{i};
                spikeSig = spikeSG.signals{i};
                SG.signals{end+1} = [baseSig(:,1) baseSig(:,2)+spikeSig(:,2)];
            end
        end
        
        function SG = gen_time_restricted_hunting(SG, number)
            huntingSG = copy(SG);
            huntingSG.type = 'hunting';
            huntingSG = gen_signals(huntingSG, number);
            lt = SG.time_range.lower;
            ut = SG.time_range.upper;
            ld = SG.min_duration.lower;
            ud = SG.min_duration.upper;
            
            durations = ld + (ud-ld)*rand(1,number);
            lrange = lt + (ut - durations).*rand(1,number);
            urange = lrange + durations;
            for i=1:number
                sig = huntingSG.signals{i};
                sig(sig(:,1)<=lrange(i) | sig(:,1) >=urange(i), 2) = mean(sig(:,2));
                SG.signals{end+1} = sig;
            end
        end
        
        function SG = gen_hunting_along_SG(SG, number)
            baseSG = copy(SG);
            baseSG.type = SG.baseName;
            baseSG = baseSG.gen_signals(number);
            huntingSG = gen_time_restricted_hunting(SG, number);
            for i=1:number
                baseSig = baseSG.signals{i};
                huntingSig = huntingSG.signals{i};
                SG.signals{end+1} = [baseSig(:,1) baseSig(:,2)+huntingSig(:,2)];
            end
        end
        
        function SG = gen_hunting_after_edge(SG, number)
            edgeSG = copy(SG);
            edgeSG.type = SG.edgeName;
            edgeSG = gen_signals(edgeSG, number);
            huntingSG = gen_sinusoid(SG,'plain', number);
            for i=1:number
                edgeSig = edgeSG.signals{i};
                huntingSig = huntingSG.signals{i};
                z = find(abs(diff(edgeSig(:,2)))>0);
                edgeSig(z:end,2) = edgeSig(z:end,2) + huntingSig(z:end,2);
                SG.signals{end+1} = edgeSig;
            end
        end
        
        function SG = gen_shoot( SG,type, number)
            if (strcmp(type,'over'))
                edgeSG = copy(SG);
                edgeSG = gen_edge(edgeSG, 'rising', number);
                mult = 1;
            else
                edgeSG = copy(SG);
                edgeSG = gen_edge(edgeSG, 'falling',number);
                mult = -1;
            end
            
            for i=1:number
                u = edgeSG.signals{i};
                uSig = u(:,2);
                uTimes = u(:,1);
                edgeTimeIndex = find(abs(diff(uSig))>0, 1, 'first');
                edgeTime = uTimes(edgeTimeIndex);
                stepValue = mult*uSig(edgeTimeIndex+1,1);
                
                if isfield(SG,'gain')
                    gain_upper = SG.gain.upper*stepValue;
                    gain_lower = SG.gain.lower*stepValue;
                else
                    gain_upper = stepValue;
                    gain_lower = stepValue;
                end
                
                if isfield(SG,'zeta')
                    transient_cfg = struct(...
                        'stepSize', SG.step_size, ...
                        'timeHorizon', SG.time_horizon, ...
                        'name', 'underdamped', ...
                        'gain', struct('lower', gain_lower, 'upper', gain_upper), ...
                        'step_time', struct('lower', edgeTime, 'upper', edgeTime),...
                        'step_amplitude', struct('lower', 1, 'upper', 1), ...
                        'omega', struct('lower', SG.omega.lower, ...
                        'upper', SG.omega.upper),...
                        'zeta', struct('lower', SG.zeta.lower, ...
                        'upper', SG.zeta.upper));
                else
                    transient_cfg = struct(...
                        'stepSize', SG.step_size, ...
                        'timeHorizon', SG.time_horizon, ...
                        'name', 'underdamped', ...
                        'gain', struct('lower', gain_lower, 'upper', gain_upper), ...
                        'step_time', struct('lower', edgeTime, 'upper', edgeTime),...
                        'step_amplitude', struct('lower', 1, 'upper', 1), ...
                        'omega', struct('lower', SG.omega.lower, ...
                        'upper', SG.omega.upper));
                end
                
                transientSG = SignalGenerator();
                transientSG.ImportSignalConfig(transient_cfg);
                transientSG.gen_signals(1);
                
                y = transientSG.signals{1};
                SG.signals{end+1} = [uTimes y(:,2) uSig];
            end
        end
        
        function SG = gen_difference(SG, number)
            % generate two random SGs with sup-norm between them at most a
            % user-specified amount
            ldiff = SG.sup_diff.lower;
            udiff = SG.sup_diff.upper;
            numCp = SG.num_control_points;
            interpType = SG.interpolationType;
            lv = SG.value_range.lower;
            uv = SG.value_range.upper;
            for i=1:number
                cpValues = lv + (uv-lv)*rand(1,numCp);
                cpTimes = linspace(0,SG.time_horizon,numCp);
                sigTimes = 0:SG.step_size:SG.time_horizon;
                if (strcmp(interpType, 'pconst'))
                    p = @(ts,xs,t) xs(min([find(ts-t>0,1)-1 length(xs)]));
                    sig1 = arrayfun(@(t) p(cpTimes, s, t), sigTimes);
                else
                    sig1 = interp1(cpTimes, cpValues, sigTimes, interpType);
                end
                incDec = -1 + 2*randi([0 1], 1, length(sigTimes));
                diffs = incDec.*(ldiff + (udiff-ldiff)*rand(1,length(sigTimes)));
                sig2 = sig1 + diffs;
                SG.signals{end+1} = [tranpose(sigTimes) tranpose(sig1) transpose(sig2)];
            end
        end
        
        function SG = gen_sticking(SG, number)
            constSG = gen_constant(SG,number);
            l = SG.value_range.lower;
            u = SG.value_range.upper;
            lr = SG.decay_rate_range.lower;
            ur = SG.decay_rate_range.upper;
            rates = lr + (ur-lr)*rand(1,number);
            vals = l + (u-l)*rand(1,number);
            
            for i=1:number
                sig1 = constSG.signals{end+1};
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
                SG.signals{end+1} = [sigTimes constSig sig2];
            end
        end
        
        function SG = gen_ramp(SG, number)
            % assumes slope != 0
            l = SG.value_range.lower;
            u = SG.value_range.upper;
            ls = SG.ramp_slope.lower;
            us = SG.ramp_slope.upper;
            lt = SG.ramp_start_time.lower;
            ut = SG.ramp_start_time.upper;
            la = SG.min_amplitude.lower;
            ua = SG.min_amplitude.upper;
            slopes = ls + (us-ls)*rand(1, number);
            startTimes = lt + (ut-lt)*rand(1,number);
            amps = la + (ua-la)*rand(1,number);
            
            lvals = l + (u-amps-l).*rand(1,number);
            uvals = lvals + amps;
            times = 0:SG.step_size:SG.time_horizon;
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
                SG.signals{end+1} = [transpose(times) transpose(sig)];
            end
        end
        
        function SG = gen_sawtooth(SG,number)
            l = SG.value_range.lower;
            u = SG.value_range.upper;
            lam = SG.min_amplitude.lower;
            uam = SG.min_amplitude.upper;
            lw = SG.width_range.lower;
            uw = SG.width_range.upper;
            lp = SG.period_range.lower;
            up = SG.period_range.upper;
            la = SG.offset_range.lower;
            ua = SG.offset_range.upper;
            amps = lam + (uam-lam)*rand(1,number);
            lvals = l + (u - uam - l).*rand(1,number);
            uvals = lvals + amps;
            times = 0:SG.step_size:SG.time_horizon;
            offsets = la + (ua-la)*rand(1,number);
            widths = lw + (uw-lw)*rand(1,number); % 0<w<1
            periods = lp + (up-lp)*rand(1,number);
            sawt = @(l,u,o,p,w,t) (l + (u-l)*sawtooth(2*pi*(t - o)/p,w));
            for i=1:number
                sig = sawt(lvals(i), uvals(i), offsets(i), ...
                    periods(i), widths(i), times);
                SG.signals{end+1} = [transpose(times) transpose(sig)];
            end
        end
        
        function SG = gen_timed_causal(SG, number)
            % alw( flag1  => ev_[0,t] flag2)
            lt = SG.time_diff.lower;
            ut = SG.time_diff.upper;
            cpArray = linspace(0, SG.time_horizon, SG.numCp);
            times = 0:SG.step_size:SG.time_horizon;
            p = @(ts,xs,t) xs(min([find(ts-t>0,1)-1 length(xs)]));
            for i=1:number
                if (rand(1,1) < 0.01)
                    % 1 % chance of producing a SG in which antecedent is never
                    % true
                    sig1cp = zeros(1,length(cpArray));
                    sig2cp = randi([0 1],1,length(cpArray));
                    sig2 = arrayfun(@(t) p(cpArray, sig2cp, t), times);
                    sig1 = arrayfun(@(t) p(cpArray, sig1cp, t), times);
                else
                    sig1cp = randi([0 1], 1, length(cpArray));
                    sig1 = arrayfun(@(t) p(cpArray, sig1cp, t), times);
                    sig2 = zeros(1,length(times));
                    endTime = SG.time_horizon - ut;
                    lastInd = find(times < endTime, 1, 'last');
                    matches = find(sig1(1,1:lastInd));
                    dt = lt + (ut-lt)*rand(1,length(matches));
                    for j=1:length(matches)
                        sig1Index = matches(j);
                        sig2Index = find(times < (times(sig1Index) + dt(j)), 1, 'last');
                        % stay high till next cp. so, first find next cp
                        sig2Last = ceil (ceil((sig2Index / length(times))*SG.numCp) * ...
                            (length(times)/SG.numCp));
                        sig2(sig2Index:sig2Last) = 1;
                    end
                    sig1(lastInd+1:end) = 0;
                end
                SG.signals{end+1} = [transpose(times) transpose(sig2(1:numel(times))) transpose(sig1(1:numel(times))) ];
            end
        end
        
        function SG = gen_pulse(SG, number)
            % assume: ut + ud < time_horizon
            lb = SG.base_value.lower;
            ub = SG.base_value.upper;
            
            la = SG.amplitude.lower;
            ua = SG.amplitude.upper;
            
            ld = SG.duration.lower;
            ud = SG.duration.upper;
            
            lt = SG.start_time.lower;
            ut = SG.start_time.upper;
            
            times = 0:SG.step_size:SG.time_horizon;
            startTimes = lt + (ut-lt)*rand(1,number);
            durations = ld + (ud-ld)*rand(1,number);
            bases = lb + (ub-lb)*rand(1,number);
            amplitudes = la + (ua-la)*rand(1,number);
            
            for i=1:number
                sig = bases(i)*ones(1,length(times));
                sig( times > startTimes(i) & ...
                    times < startTimes(i) + durations(i) ) = ...
                    bases(i) + amplitudes(i);
                SG.signals{end+1} = [transpose(times) transpose(sig)];
            end
        end
        
        function SG = gen_pulse_train(SG, number)
            % assume: ut + ud < time_horizon
            lb = SG.base_value.lower;
            ub = SG.base_value.upper;
            
            la = SG.amplitude.lower;
            ua = SG.amplitude.upper;
            
            ld = SG.duty_cycle.lower;
            ud = SG.duty_cycle.upper;
            
            lp = SG.period.lower;
            up = SG.period.upper;
            
            times = 0:SG.step_size:SG.time_horizon;
            periods = lp + (up-lp)*rand(1,number);
            dutyCycles = ld + (ud-ld)*rand(1,number);
            bases = lb + (ub-lb)*rand(1,number);
            amplitudes = la + (ua-la)*rand(1,number);
            
            unshifted_pulse = @(T,d,t) sin(((2*pi)/T)*t) > cos(pi*d);
            advance = @(T,d) T*(0.25 - 0.5*d);
            pulse = @(T,d,t) unshifted_pulse(T,d, t + advance(T,d));
            
            for i=1:number
                sig = bases(i) + amplitudes(i)*pulse(periods(i), dutyCycles(i), times);
                SG.signals{end+1} = [transpose(times) transpose(sig)];
            end
        end
        
        function SG = gen_second_order_step_response(SG,type, number)
            la = SG.step_amplitude.lower;
            ua = SG.step_amplitude.upper;
            ls = SG.step_time.lower;
            us = SG.step_time.upper;
            lw = SG.omega.lower;
            uw = SG.omega.upper;
            lg = SG.gain.lower;
            ug = SG.gain.upper;
            omegas = lw + (uw-lw)*rand(1,number);
            gains = lg + (ug-lg)*rand(1,number);
            amps = la + (ua-la)*rand(1,number);
            sts = ls + (us-ls)*rand(1,number);
            if isfield(SG,'zeta')
                zetal = SG.zeta.lower;
                zetau = SG.zeta.upper;
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
            
            times = 0:SG.step_size:SG.time_horizon;
            for i=1:number
                sys = tf([gains(i)*omegas(i)^2], [1 2*zetas(i)*omegas(i) omegas(i)^2]);
                u = amps(i)*heaviside(times-sts(i));
                SG.signals{end+1} = [transpose(times) lsim(sys,u,times)];
            end
        end
        
        function plot_signals(SG)
            plot_signals(SG.signals, SG.type);
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
        
        function disp(SG)
            
            disp('Signal Generator with Configurations:')
            disp('-------------------------------------')
            for ip = 1:numel(SG.configs)
                fprintf('%i: %s\n', ip, SG.configs{ip}.name);
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

