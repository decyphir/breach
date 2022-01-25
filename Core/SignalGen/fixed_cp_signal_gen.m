classdef fixed_cp_signal_gen < signal_gen
    % fixed_cp_signal_gen  A class derived from signal_gen to generate signals from regularly spaced control points using different interpolation methods.
    %
    % fixed_cp_signal_gen Properties
    %   cp     - number of control points
    %   method - interpolation method for each signal, accepts all methods provided by interp1
    %
    % fixed_cp_signal_gen Methods
    %   fixed_cp_signal_gen -  constructor, takes signal names, number of
    %                          control points, interpolation methods and a
    %                          p0.
    %
    %  See also signal_gen.
    
    properties
        num_cp       % number of control points for each signal
        method   % interpolation method for each signal
    end
    
    methods
        
        function this = fixed_cp_signal_gen(signals, cp, method,p0)
            
            if nargin == 1
                cp = 1;
                method = 'previous';
            elseif nargin == 2
                method = 'previous';
            end
            
            if ~iscell(signals)
                signals = {signals};
            end
            
            if isscalar(cp)
                cp = repmat(cp,1, numel(signals));
            end
            
            if ~iscell(method)
                method = {method};
            end
            
            
            this.signals = signals;
            this.num_cp = cp;
            this.params = {};
            
            if nargin<=2
                method = 'previous';
            end
            
            this.method = method;
            for ku = 1:numel(signals)
                for k = 1:cp(ku)
                    this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
                end
            end
            if nargin==4
                this.p0 = p0;
            else
                this.p0= zeros(numel(this.params),1);
            end
            
            this.params_domain = repmat(BreachDomain(), 1, numel(this.params));
            this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));
            
        end
        
        function [X, time] = computeSignals(this,p,time) % compute the signals
            if size(p,1) ==1
                p = p';
            end
            
            X = zeros(numel(this.num_cp),numel(time));
            pts_x = p;
            for i_cp = 1:numel(this.num_cp)
                cp_values = pts_x(1:this.num_cp(i_cp));
                pts_x = pts_x(this.num_cp(i_cp)+1:end);
                if ischar(this.method)
                    meth = this.method;
                elseif numel(this.method)==1
                    meth = this.method{1};
                else
                    meth = this.method{i_cp};
                end
                switch meth
                    case 'previous'
                        t_cp = linspace(time(1), time(end), this.num_cp(i_cp)+1)';
                        if numel(t_cp)==2
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            x = interp1(t_cp(1:end-1), cp_values, time', meth, 'extrap');
                        end
                        X(i_cp,:) = x';
                        
                    otherwise
                        t_cp = linspace(time(1), time(end), this.num_cp(i_cp))';
                        if numel(t_cp) == 1
                            x = cp_values(1)*ones(numel(time),1);
                        else
                            x = interp1(t_cp, cp_values, time', meth, 'extrap');
                        end
                        X(i_cp,:) = x';
                end
            end
        end
        function type = getType(this)
            type = 'unistep';
        end
        
        function args = getSignalGenArgs(this)
            args = {'num_cp','method'};
        end
        
        function [t_cp, cp_values, cp_low, cp_up, dom] = get_cp(this, signal, time)
            
            i_sig = find(strcmp(signal,this.signals),1);
            i_p = sum(this.num_cp(1:i_sig-1))+1;
            idx_cp = i_p:i_p+this.num_cp(i_sig)-1;
            
            cp_values = this.p0(idx_cp);
            if ischar(this.method)
                meth = this.method;
            elseif numel(this.method)==1
                meth = this.method{1};
            else
                meth = this.method{i_cp};
            end
            switch meth
                case 'previous'
                    t_cp = linspace(time(1), time(end), this.num_cp(i_sig)+1)';
                    t_cp = t_cp(1:end-1);
                otherwise
                    t_cp = linspace(time(1), time(end), this.num_cp(i_sig))';
            end
            
            cp_low = cp_values;
            cp_up = cp_values;
            
            
            for ic = 1:numel(idx_cp)
                dom(ic) = this.params_domain(idx_cp(ic));
                if ~isempty(dom(ic).domain)
                    cp_low(ic) = dom(ic).domain(1);
                    cp_up(ic) = dom(ic).domain(2);
                end
            end
            
        end
        
        function plot(this, signal, time)
            plot@signal_gen(this,signal, time);
            
            % plot control points
            [t_cp, x_cp] = this.get_cp(signal, time);
            hold on;
            plot(t_cp,x_cp,'or','MarkerSize', 6, 'MarkerFaceColor', [1 0 0]);
            legend({signal,'control points (cp)'}, 'Interpreter', 'None');
            
        end
        
        function plot_enveloppe(this, signal, time, varargin)
            
            
            opt.max_time = 3;
            opt.max_obj_eval = inf;
            opt.constraints = {}; 
            opt = varargin2struct_breach(opt, varargin{:});
                        
            isig= find(strcmp(signal,this.signals),1);
            
            if iscell(this.method)
                method = this.method{isig};
            else
                method = this.method;
            end
            
            if ~isempty(opt.constraints)||~ismember(method, {'previous', 'linear'});
                plot_enveloppe@signal_gen(this,signal, time, varargin{:});
            else
                
                [t_cp, x_cp, x_low, x_up, dom] = this.get_cp(signal, time);
                
                
                switch (method)
                    case 'previous'
                        for icp = 1:numel(t_cp)-1
                            %switch (dom(ic).type)
                            %    case 'double'
                            X = [t_cp(icp) t_cp(icp) t_cp(icp+1) t_cp(icp+1)];
                            Y = [x_low(icp) x_up(icp) x_up(icp) x_low(icp)];
                            patch(X, Y, 'k', 'FaceAlpha', 0.1, 'EdgeAlpha', .1);
                        end
                        X = [t_cp(end) t_cp(end) time(end) time(end)];
                        Y = [x_low(end) x_up(end) x_up(end) x_low(end)];
                        patch(X, Y, 'k', 'FaceAlpha', 0.1, 'EdgeAlpha',.1);
                        
                        
                    otherwise
                        
                        for icp = 1:numel(t_cp)-1
                            %switch (dom(ic).type)
                            %    case 'double'
                            X = [t_cp(icp) t_cp(icp) t_cp(icp+1) t_cp(icp+1)];
                            Y = [x_low(icp) x_up(icp) x_up(icp+1) x_low(icp+1)];
                            patch(X, Y, 'k', 'FaceAlpha', 0.1, 'EdgeAlpha', .1);
                        end
                        X = [t_cp(end) t_cp(end) time(end) time(end)];
                        Y = [x_low(end) x_up(end) x_up(end) x_low(end)];
                        patch(X, Y, 'k', 'FaceAlpha', 0.1, 'EdgeAlpha', .1);
                end
                
            end
            
        end
    end
end
