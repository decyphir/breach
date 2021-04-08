classdef testron_signal_gen < signal_gen
    properties
        args
        signal_type % some parameter for the signal generator - this won't be visible from Breach API
    end
    methods
        
        function this = testron_signal_gen(signals, signal_type, args)
            % Args can be the following:
            
            % signal_type = discrete_constant
            % args = [min_values, max_values, n_intervals]
            
            % signal_type = continuous_constant
            % args = [min_values, max_values, n_intervals]
            
            this.signals = signals; % Signal names
            this.signal_type = signal_type; % Signal type
            this.args = args;
            this.params = {};
            this.p0 = [];
            
            for ku = 1:numel(signals)
                if strcmp(signal_type{ku},'discrete_constant')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'discrete_constant_fixed_start')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    start_value = this_arg(4);
                    this.params = {this.params{:} [signals{ku} '_u0']};
                    this.p0(end+1) = start_value;
                    for k = 2:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'discrete_constant_fixed_start_interval')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'continuous_constant')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'continuous_constant_fixed_start_interval')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'continuous_pulse')
                    this_arg = this.args{ku};
                    min_amp = this_arg(1);
                    max_amp = this_arg(2);
                    min_period = this_arg(3);
                    max_period = this_arg(4);
                    base_val = (min_amp + max_amp)/2;
                    
                    % Add parameter names
                    this.params = { this.params{:} [this.signals{ku} '_base_value'] ...
                        [this.signals{ku} '_pulse_period']...
                        [this.signals{ku} '_pulse_width']...
                        [this.signals{ku} '_pulse_amp']...
                        [this.signals{ku} '_pulse_delay']};
                    
                    % Set standard parameter values
                    this.p0(end+1) = 0; % Base value
                    this.p0(end+1) = min_period; % Pulse period (will be set by user)
                    this.p0(end+1) = 0.5; % Pulse width
                    this.p0(end+1) = min_amp; % Pulse amp (will be set by user)
                    if min_amp == max_amp
                        this.p0(end+1) = 5;
                    else
                        this.p0(end+1) = 0; % Pulse delay
                    end
                    
                elseif strcmp(signal_type{ku},'discrete_pulse')
                    this_arg = this.args{ku};
                    min_amp = this_arg(1);
                    max_amp = this_arg(2);
                    min_period = this_arg(3);
                    max_period = this_arg(4);
                    base_val = (min_amp + max_amp)/2;
                    
                    % Add parameter names
                    this.params = { this.params{:} [this.signals{ku} '_base_value'] ...
                        [this.signals{ku} '_pulse_period']...
                        [this.signals{ku} '_pulse_width']...
                        [this.signals{ku} '_pulse_amp']...
                        [this.signals{ku} '_pulse_delay']};
                    
                    % Set standard parameter values
                    this.p0(end+1) = 0; % Base value
                    this.p0(end+1) = min_period; % Pulse period (will be set by user)
                    this.p0(end+1) = 0.5; % Pulse width
                    this.p0(end+1) = min_amp; % Pulse amp (will be set by user)
                    if min_amp == max_amp
                        this.p0(end+1) = 5;
                    else
                        this.p0(end+1) = 0; % Pulse delay
                    end
                    
                elseif strcmp(signal_type{ku},'continuous_cp')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'discrete_cp')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'continuous_cp_fixed_start')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    start_value = this_arg(4);
                    this.params = {this.params{:} [signals{ku} '_u0']};
                    this.p0(end+1) = start_value;
                    for k = 2:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'discrete_cp_fixed_start')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    start_value = this_arg(4);
                    this.params = {this.params{:} [signals{ku} '_u0']};
                    this.p0(end+1) = start_value;
                    for k = 2:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'cp_set_times')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(1);
                    cp_times = this_arg(2:2+n_intervals-1);
                    cp_values = this_arg(2+n_intervals:end);
                    for k = 1:n_intervals
                        % Add the value
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = cp_values(k);
                        
                        % Add the time
                        this.params= {this.params{:} [signals{ku} '_t' num2str(k-1)]};
                        this.p0(end+1) = cp_times(k);
                    end
                    
                elseif strcmp(signal_type{ku},'No_GearEven')
                    this_arg = this.args{ku};
                    this.params= {this.params{:} [signals{ku} '_u0']};
                    this.p0(end+1) = 0;
                    
                elseif strcmp(signal_type{ku},'continuous_cp_fixed_start_period')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    start_value = this_arg(4);
                    start_period_time = this_arg(5);
                    this.params = {this.params{:} [signals{ku} '_u0']};
                    this.p0(end+1) = start_value;
                    for k = 2:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'discrete_cp_fixed_start_period')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    start_value = this_arg(4);
                    start_period_time = this_arg(5);
                    this.params = {this.params{:} [signals{ku} '_u0']};
                    this.p0(end+1) = start_value;
                    for k = 2:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'continuous_cp_fixed_start_interval')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'discrete_cp_fixed_start_interval')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(3);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                elseif strcmp(signal_type{ku},'discrete_enumeration')
                    % Previously No_GearEven
                    this_arg = this.args{ku};
                    n_intervals = this_arg(end);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                
                elseif strcmp(signal_type{ku},'discrete_enumeration_fixed_start_period')
                    this_arg = this.args{ku};
                    n_intervals = this_arg(end);
                    for k = 1:n_intervals
                        this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                        this.p0(end+1) = 0;
                    end
                    
                else
                    disp(['Warning! Cannot handle signal type ' signal_type]);
                end
                
            end
            
        end
        
        function X = computeSignals(this,p,time) % compute the signals
            if size(p,1) ==1
                p = p';
            end
            
            X = zeros(numel(this.signals),numel(time));
            pts_x = p;
            
            % Debugging
            %param_counter = evalin('base','param_counter;');
            %assignin('base',['param' num2str(param_counter)],pts_x);
            %evalin('base','param_counter = param_counter + 1;');
            %interesting_idx = find(strcmp(this.signals, 'IsgModReq__IsgModReq'));
            interesting_idx = -100;
            
            % Loop over each signal
            for i_ni = 1:numel(this.signals)
                
                switch this.signal_type{i_ni}
                    
                    case 'discrete_constant'
                        this_arg = this.args{i_ni};
                        max_value = this_arg(2);
                         min_value = this_arg(1);
                       n_intervals = this_arg(3);
                        interval_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        interval_values = floor(interval_values);
                        t_cp = linspace(time(1), time(end), n_intervals+1)';
                        if numel(t_cp)==1
                            x = interval_values(1)*ones(numel(time),1);
                        else
                            x = zeros(size(time));
                            for tmp = 1:length(time)-1
                                interval_index = find(t_cp > time(tmp),1) - 1;
                                x(tmp) = interval_values(interval_index);
                            end
                            x(end) = interval_values(end);
                            x = x';
                        end
                        x = min(x,max_value);
                        x = max(x, min_value);
                        X(i_ni,:) = x';
                        
                    case 'discrete_constant_fixed_start'
                        this_arg = this.args{i_ni};
                        min_value = this_arg(1);
                        max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        start_value = this_arg(4);
                        interval_values = pts_x(1:n_intervals);
                        interval_values = [0; interval_values];
                        pts_x = pts_x(n_intervals+2:end);
                        interval_values = floor(interval_values);
                        t_cp = linspace(time(1), time(end), n_intervals+2)';
                        if numel(t_cp)==1
                            x = interval_values(1)*ones(numel(time),1);
                        else
                            x = zeros(size(time));
                            for tmp = 1:length(time)-1
                                interval_index = find(t_cp > time(tmp),1) - 1;
                                x(tmp) = interval_values(interval_index);
                            end
                            x(end) = interval_values(end);
                            x = x';
                        end
                        x = min(x,max_value);
                        x = max(x, min_value);
                        X(i_ni,:) = x';
                        
                    case 'discrete_constant_fixed_start_interval'
                        this_arg = this.args{i_ni};
                        max_value = this_arg(2);
                        min_value = this_arg(1);
                        n_intervals = this_arg(3);
                        interval_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        interval_values = floor(interval_values);
                        t_cp = linspace(time(1), time(end), n_intervals+1)';
                        if numel(t_cp)==1
                            x = interval_values(1)*ones(numel(time),1);
                        else
                            x = zeros(size(time));
                            for tmp = 1:length(time)-1
                                interval_index = find(t_cp > time(tmp),1) - 1;
                                x(tmp) = interval_values(interval_index);
                            end
                            x(end) = interval_values(end);
                            x = x';
                        end
                        x = min(x,max_value);
                        x = max(x, min_value);
                        X(i_ni,:) = x';
                        
                    case 'continuous_constant'
                        this_arg = this.args{i_ni};
                        %min_value = this_arg(1);
                        %max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        interval_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(time(1), time(end), n_intervals+1)';
                        if numel(t_cp)==1
                            x = interval_values(1)*ones(numel(time),1);
                        else
                            x = zeros(size(time));
                            for tmp = 1:length(time)-1
                                interval_index = find(t_cp > time(tmp),1) - 1;
                                x(tmp) = interval_values(interval_index);
                            end
                            x(end) = interval_values(end);
                            x = x';
                        end
                        X(i_ni,:) = x';
                        
                    case 'continuous_constant_fixed_start_interval'
                        this_arg = this.args{i_ni};
                        %min_value = this_arg(1);
                        %max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        interval_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(time(1), time(end), n_intervals+1)';
                        if numel(t_cp)==1
                            x = interval_values(1)*ones(numel(time),1);
                        else
                            x = zeros(size(time));
                            for tmp = 1:length(time)-1
                                interval_index = find(t_cp > time(tmp),1) - 1;
                                x(tmp) = interval_values(interval_index);
                            end
                            x(end) = interval_values(end);
                            x = x';
                        end
                        X(i_ni,:) = x';
                        
                    case 'continuous_pulse'
                        unshifted_pulse = @(T,d,t) sin(((2*pi)/T)*(t+eps)) > cos(pi*d);
                        advance = @(T,d) T*(0.25 - 0.5*d);
                        pulse = @(T,d,t) unshifted_pulse(T,d, t + advance(T,d));
                        
                        base = pts_x(1);
                        period = pts_x(2);
                        duty = pts_x(3);
                        amp = pts_x(4);
                        delay = pts_x(5);
                        pts_x = pts_x(6:end);
                        
                        if delay ==0
                            X(i_ni,:) = base + amp*pulse(period, duty, time);
                        elseif delay>0
                            X(i_ni, :) = amp*ones(size(time));
                        else
                            
                        end
                        
                    case 'discrete_pulse'
                        unshifted_pulse = @(T,d,t) sin(((2*pi)/T)*(t+eps)) > cos(pi*d);
                        advance = @(T,d) T*(0.25 - 0.5*d);
                        pulse = @(T,d,t) unshifted_pulse(T,d, t + advance(T,d));
                        
                        base = pts_x(1);
                        period = pts_x(2);
                        duty = pts_x(3);
                        amp = pts_x(4);
                        delay = pts_x(5);
                        pts_x = pts_x(6:end);
                        
                        if delay ==0
                            X(i_ni,:) = floor(base + amp*pulse(period, duty, time));
                        elseif delay>0
                            X(i_ni, :) = amp*ones(size(time));
                        else
                            
                        end
                        
                    case 'continuous_cp'
                        this_arg = this.args{i_ni};
                        min_value = this_arg(1);
                        max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        cp_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(time(1), time(end), n_intervals)';
                        if numel(t_cp)==1
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            %x = interp1(t_cp, cp_values, time', 'linear', 'extrap');
                            x = interp1(t_cp, cp_values, time', 'pchip', 'extrap');
                        end
                        x = min(x,max_value);
                        x = max(x, min_value);
                        X(i_ni,:) = x';
                    
                    case 'discrete_cp'
                        
                        this_arg = this.args{i_ni};
                        min_value = this_arg(1);
                        max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        cp_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(time(1), time(end), n_intervals)';
%                         if i_ni == interesting_idx
%                             %disp('Warning! Changing the input generated for IsgModReq__IsgModReq');
%                             min_value = this_arg(1);
%                             max_value = this_arg(2) + 1;
%                             r = (max_value-min_value).*rand(n_intervals,1) + min_value;
%                             cp_values = floor(r);
%                         end
                        if numel(t_cp)==1
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            %x = interp1(t_cp, cp_values, time', 'linear', 'extrap');
                            x = interp1(t_cp, cp_values, time', 'pchip', 'extrap');
                        end
                        x = floor(x);
                        x = min(x,max_value);
                        x = max(x, min_value);
                        X(i_ni,:) = x';
                        
                    case 'continuous_cp_fixed_start'
                        this_arg = this.args{i_ni};
                        min_value = this_arg(1);
                        max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        start_value = this_arg(4);
                        cp_values = pts_x(1:n_intervals);
                        cp_values(1) = start_value;
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(time(1), time(end), n_intervals)';
                        if numel(t_cp)==1
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            x = interp1(t_cp, cp_values, time', 'pchip', 'extrap');
                        end
                        x = min(x,max_value);
                        x = max(x, min_value);
                        X(i_ni,:) = x';
                        
                        
                    case 'discrete_cp_fixed_start'
                        
                        this_arg = this.args{i_ni};
                        min_value = this_arg(1);
                        max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        start_value = this_arg(4);
                        cp_values = pts_x(1:n_intervals);
                        cp_values(1) = start_value;
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(time(1), time(end), n_intervals)';
                        if i_ni == interesting_idx
                            %disp('Warning! Changing the input generated for IsgModReq__IsgModReq');
                            min_value = this_arg(1);
                            max_value = this_arg(2) + 1;
                            r = (max_value-min_value).*rand(n_intervals,1) + min_value;
                            cp_values = floor(r);
                        end
                        if numel(t_cp)==1
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            %x = interp1(t_cp, cp_values, time', 'linear', 'extrap');
                            x = interp1(t_cp, cp_values, time', 'pchip', 'extrap');
                        end
                        x = floor(x);
                        x = min(x,max_value);
                        x = max(x, min_value);
                        X(i_ni,:) = x';
                        
                    case 'cp_set_times'
                        this_arg = this.args{i_ni};
                        n_intervals = this_arg(1);
                        cp_times = this_arg(2:2+n_intervals-1);
                        cp_values = this_arg(2+n_intervals:end);
                        pts_x = pts_x(1 + 2*n_intervals:end);
                        x = zeros(numel(time),1);
                        for k = 1:n_intervals
                            x(time >= cp_times(k)) = cp_values(k);  
                        end
                        X(i_ni,:) = x';
                        
                    case 'No_GearEven'
                        this_arg = this.args{i_ni};
                        %min_value = this_arg(1);
                        %max_value = this_arg(2);
                        n_intervals = 1;
                        interval_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        interval_values = floor(interval_values);
                        try
                            interval_values = this_arg(interval_values);
                        catch
                            % The value of interval_values was EXACTLY 6
                            interval_values = this_arg(end);
                        end
                        t_cp = linspace(time(1), time(end), n_intervals+1)';
                        if numel(t_cp)==1
                            x = interval_values(1)*ones(numel(time),1);
                        else
                            x = zeros(size(time));
                            for tmp = 1:length(time)-1
                                interval_index = find(t_cp > time(tmp),1) - 1;
                                x(tmp) = interval_values(interval_index);
                            end
                            x(end) = interval_values(end);
                            x = x';
                        end
                        X(i_ni,:) = x';
                        
                    case 'continuous_cp_fixed_start_period'
                        % Fixed start, but at a DELAYED start time!
                        this_arg = this.args{i_ni};
                        min_value = this_arg(1);
                        max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        start_value = this_arg(4);
                        start_period_time = this_arg(5);
                        cp_values = pts_x(1:n_intervals);
                        cp_values(1) = start_value;
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(start_period_time, time(end), n_intervals)';
                        if numel(t_cp)==1
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            x = interp1(t_cp, cp_values, time', 'pchip', 'extrap');
                        end
                        x = min(x,max_value);
                        x = max(x, min_value);
                        
                        start_period_end_index = find(time > start_period_time, 1);
                        x(1:start_period_end_index) = start_value;
                        X(i_ni,:) = x';
                        
                    case 'discrete_cp_fixed_start_period'
                        % Fixed start, but at a DELAYED start time!
                        this_arg = this.args{i_ni};
                        min_value = this_arg(1);
                        max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        start_value = this_arg(4);
                        start_period_time = this_arg(5);
                        cp_values = pts_x(1:n_intervals);
                        cp_values(1) = start_value;
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(start_period_time, time(end), n_intervals)';
                        if numel(t_cp)==1
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            %x = interp1(t_cp, cp_values, time', 'linear', 'extrap');
                            x = interp1(t_cp, cp_values, time', 'pchip', 'extrap');
                        end
                        x = floor(x);
                        
                        x = min(x,max_value);
                        x = max(x, min_value);
                        
                        start_period_end_index = find(time > start_period_time, 1);
                        x(1:start_period_end_index) = start_value;
                        X(i_ni,:) = x';
                        
                    case 'continuous_cp_fixed_start_interval'
                        % Fixed start, but from an INTERVAL instead of a
                        % single value!
                        this_arg = this.args{i_ni};
                        min_value = this_arg(1);
                        max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        cp_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(time(1), time(end), n_intervals)';
                        if numel(t_cp)==1
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            x = interp1(t_cp, cp_values, time', 'pchip', 'extrap');
                        end
                        x = min(x,max_value);
                        x = max(x, min_value);
                        X(i_ni,:) = x';
                        
                    case 'discrete_cp_fixed_start_interval'
                        % Fixed start, but from an INTERVAL instead of a
                        % single value!
                        this_arg = this.args{i_ni};
                        min_value = this_arg(1);
                        max_value = this_arg(2);
                        n_intervals = this_arg(3);
                        cp_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        t_cp = linspace(time(1), time(end), n_intervals)';
%                         if i_ni == interesting_idx
%                             %disp('Warning! Changing the input generated for IsgModReq__IsgModReq');
%                             min_value = this_arg(1);
%                             max_value = this_arg(2) + 1;
%                             r = (max_value-min_value).*rand(n_intervals,1) + min_value;
%                             cp_values = floor(r);
%                         end
                        if numel(t_cp)==1
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            %x = interp1(t_cp, cp_values, time', 'linear', 'extrap');
                            x = interp1(t_cp, cp_values, time', 'pchip', 'extrap');
                        end
                        x = floor(x);
                        x = min(x,max_value);
                        x = max(x, min_value);
                        X(i_ni,:) = x';
                        
                    case 'discrete_enumeration'
                        this_arg = this.args{i_ni};
                        %min_value = this_arg(1);
                        %max_value = this_arg(2);
                        n_intervals = this_arg(end);
                        interval_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        interval_values = floor(interval_values);
                        for k = 1:n_intervals
                            if interval_values(k) == length(this_arg)
                                interval_values(k) = interval_values(k) - 1;
                            end
                            interval_values(k) = this_arg(interval_values(k));
                        end
                        t_cp = linspace(time(1), time(end), n_intervals+1)';
                        if numel(t_cp)==1
                            x = interval_values(1)*ones(numel(time),1);
                        else
                            x = zeros(size(time));
                            for tmp = 1:length(time)-1
                                interval_index = find(t_cp > time(tmp),1) - 1;
                                x(tmp) = interval_values(interval_index);
                            end
                            x(end) = interval_values(end);
                            x = x';
                        end
                        X(i_ni,:) = x';
                        
                    case 'discrete_enumeration_fixed_start_period'
                        this_arg = this.args{i_ni};
                        start_value = this_arg(end-2);
                        start_period_time = this_arg(end-1);
                        n_intervals = this_arg(end);
                        interval_values = pts_x(1:n_intervals);
                        pts_x = pts_x(n_intervals+1:end);
                        interval_values = floor(interval_values);
                        for k = 1:n_intervals
                            if interval_values(k) == length(this_arg)
                                interval_values(k) = interval_values(k) - 1;
                            end
                            interval_values(k) = this_arg(interval_values(k));
                        end
                        t_cp = linspace(time(1), time(end), n_intervals+1)';
                        if numel(t_cp)==1
                            x = interval_values(1)*ones(numel(time),1);
                        else
                            x = zeros(size(time));
                            for tmp = 1:length(time)-1
                                interval_index = find(t_cp > time(tmp),1) - 1;
                                x(tmp) = interval_values(interval_index);
                            end
                            x(end) = interval_values(end);
                            x = x';
                        end
                        
                        start_period_end_index = find(time > start_period_time, 1);
                        x(1:start_period_end_index) = start_value;
                        X(i_ni,:) = x';
                        
                        
                    otherwise
                        error('Error! Signal type is not recognized');
                end
            end
            % TODO!
            5;
            
            % If we are generating for the CMA_CIDD model, apply input
            % constraints!
            if this.modelIsCMACIDD()
                X = this.applyCMACIDDConstraints(X, time);
            end
            
            
            
        end
        
        function type = getType(this)
            type = 'unistep';
        end
        
        function isCMACIDD = modelIsCMACIDD(this)
            % Check that all signals that should be constrained ACTUALLY
            % EXIST in the system!
            HvTestIdx = getSignalIndex(this, 'HvTestSts__HvTestSts');
            LvTestIdx = getSignalIndex(this, 'LvTestSts__LvTestSts');
            IsgModIdx = getSignalIndex(this, 'IsgModReq__IsgModReq');
            CloseClutchIdx = getSignalIndex(this, 'B_CloseClutchEven');
            SpeedControlIdx = this.getSignalIndex('B_EngineSpeedControl');
            EngineSetPointIdx = this.getSignalIndex('W_Engine_SetPoint');
            GearEvenIdx = this.getSignalIndex('No_GearEven');
            
            if ~isempty(HvTestIdx) && ...
                    ~isempty(LvTestIdx) && ...
                    ~isempty(IsgModIdx) && ...
                    ~isempty(CloseClutchIdx) && ...
                    ~isempty(SpeedControlIdx) && ...
                    ~isempty(EngineSetPointIdx) && ...
                    ~isempty(GearEvenIdx)
                isCMACIDD = 1;
            else
                isCMACIDD = 0;
            end
        end
        
        function X = applyCMACIDDConstraints(this, X, time)
            % Load all the signal indexes we need
            HvTestIdx = this.getSignalIndex('HvTestSts__HvTestSts');
            LvTestIdx = this.getSignalIndex('LvTestSts__LvTestSts');
            IsgModIdx = this.getSignalIndex('IsgModReq__IsgModReq');
            CloseClutchIdx = this.getSignalIndex('B_CloseClutchEven');
            SpeedControlIdx = this.getSignalIndex('B_EngineSpeedControl');
            EngineSetPointIdx = this.getSignalIndex('W_Engine_SetPoint');
            GearEvenIdx = this.getSignalIndex('No_GearEven');
            
            % IsgSpdReq__IsgSpdReq
            % Depends on ModStsIsg__ModStsIsg (not an INPUT)
            % =======
            % if ModStsIsg != 3
            %   IsgSpdReq = 0
            
            
            % HvTestSts__HvTestSts AND LvTestSts__LvTestSts
            % =======
            % if time > 0.5
            %   HvTestSts = 2
            %   LvTestSts = 2
            X(HvTestIdx, time > 0.5) = 2;
            X(LvTestIdx, time > 0.5) = 2;
            
            % IsgTqReq__IsgTqReq
            % Depends on ModStsIsg__ModStsIsg (not an INPUT)
           
            % W_Engine_SetPoint
            % =======
            % if B_EngineSpeedControl == 0
            %   W_Engine_SetPoint = 0
            % else
            %   W_Engine_SetPoint >= 100
            X(EngineSetPointIdx, X(SpeedControlIdx,:) == 0) = 0;
            X(EngineSetPointIdx, X(SpeedControlIdx,:) ~= 0) = ...
                min(X(EngineSetPointIdx, X(SpeedControlIdx,:) ~= 0), 100);
            
            % B_CloseClutchEven
            % =======
            % if IsgModReq == 2 or IsgModReq == 4
            %   B_CloseClutchEven = 1
            % else
            %   B_CloseClutchEven = 0
            X(CloseClutchIdx, X(IsgModIdx,:) == 2 | ...
                X(IsgModIdx,:) == 4) = 1;
            X(CloseClutchIdx, X(IsgModIdx,:) ~= 2 & ...
                X(IsgModIdx,:) ~= 4) = 0;
            
            % No_GearEven
            % =======
            % if B_CloseClutchEven == 1 or IsgModReq != 2
            %   No_GearEven = 0
            % elseif IsgModReq == 2
            %   Do nothing, No_GearEven is chosen correctly
            % else
            %   No_GearEven = 0
            X(GearEvenIdx, X(CloseClutchIdx,:) == 1 | ...
                X(IsgModIdx,:) ~=2) = 0;
        end
        
        function isV331IEM = modelIsV331IEM(this)
            % Check that all signals that should be constrained ACTUALLY
            % EXIST in the system!
            HvTestIdx = getSignalIndex(this, 'HvTestSts__HvTestSts');
            LvTestIdx = getSignalIndex(this, 'LvTestSts__LvTestSts');
            CloseClutchIdx = getSignalIndex(this, 'B_CloseClutchEven');
            SpeedControlIdx = this.getSignalIndex('B_EngineSpeedControl');
            EngineSetPointIdx = this.getSignalIndex('W_Engine_SetPoint');
            GearEvenIdx = this.getSignalIndex('No_GearEven');
            
            if ~isempty(HvTestIdx) && ...
                    ~isempty(LvTestIdx) && ...
                    ~isempty(CloseClutchIdx) && ...
                    ~isempty(SpeedControlIdx) && ...
                    ~isempty(EngineSetPointIdx) && ...
                    ~isempty(GearEvenIdx)
                isV331IEM = 1;
            else
                isV331IEM = 0;
            end
        end
        
        function X = applyV331IEMConstraints(this, X, time)
            % This is the same as applyCMACIDDConstraints, but without
            % using IsgModReq (IsgModIdx), since it does not exist in
            % V331_IEM
            % Load all the signal indexes we need
            HvTestIdx = this.getSignalIndex('HvTestSts__HvTestSts');
            LvTestIdx = this.getSignalIndex('LvTestSts__LvTestSts');
            SpeedControlIdx = this.getSignalIndex('B_EngineSpeedControl');
            EngineSetPointIdx = this.getSignalIndex('W_Engine_SetPoint');
            
            % IsgSpdReq__IsgSpdReq
            % Depends on ModStsIsg__ModStsIsg (not an INPUT)
            % =======
            % if ModStsIsg != 3
            %   IsgSpdReq = 0
            
            
            % HvTestSts__HvTestSts AND LvTestSts__LvTestSts
            % =======
            % if time > 0.5
            %   HvTestSts = 2
            %   LvTestSts = 2
            X(HvTestIdx, time > 0.5) = 2;
            X(LvTestIdx, time > 0.5) = 2;
            
            % IsgTqReq__IsgTqReq
            % Depends on ModStsIsg__ModStsIsg (not an INPUT)
           
            % W_Engine_SetPoint
            % =======
            % if B_EngineSpeedControl == 0
            %   W_Engine_SetPoint = 0
            % else
            %   W_Engine_SetPoint >= 100
            X(EngineSetPointIdx, X(SpeedControlIdx,:) == 0) = 0;
            X(EngineSetPointIdx, X(SpeedControlIdx,:) ~= 0) = ...
                min(X(EngineSetPointIdx, X(SpeedControlIdx,:) ~= 0), 100);
            
            % B_CloseClutchEven
            % =======
            % if IsgModReq == 2 or IsgModReq == 4
            %   B_CloseClutchEven = 1
            % else
            %   B_CloseClutchEven = 0
%             X(CloseClutchIdx, X(IsgModIdx,:) == 2 | ...
%                 X(IsgModIdx,:) == 4) = 1;
%             X(CloseClutchIdx, X(IsgModIdx,:) ~= 2 & ...
%                 X(IsgModIdx,:) ~= 4) = 0;
            
            % No_GearEven
            % =======
            % if B_CloseClutchEven == 1 or IsgModReq != 2
            %   No_GearEven = 0
            % elseif IsgModReq == 2
            %   Do nothing, No_GearEven is chosen correctly
            % else
            %   No_GearEven = 0
%             X(GearEvenIdx, X(CloseClutchIdx,:) == 1 | ...
%                 X(IsgModIdx,:) ~=2) = 0;
        end
        
        function idx = getSignalIndex(this, signal)
            IndexC = strcmp(this.signals, signal);
            idx = find(IndexC > 0);
        end
    end
end