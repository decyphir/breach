classdef var_step_signal_gen < signal_gen
    % var_step_signal_gen  A class derived from signal_gen to generate signals from control points
    %                    with variable time steps, using different
    %                    interpolation methods. Difference with var_cp is
    %                    that each signal share the same time steps
    %
    % var_cp_signal_gen Properties
    %   cp     - number of control points
    %   method - interpolation method for each signal, accepts all methods provided by interp1
    %
    % var_cp_signal_gen Methods
    %   var_step_signal_gen -  constructor, takes signal names, number of
    %                          control points, interpolation methods and a
    %                          p0.
    %  See also signal_gen.
    
    % TOFIX can't combine two var_step_signal_gen, cause conflicts with
    %       dt_u0, dt_u1 names...
    properties
        cp       % number of control points for each signal
        method   % interpolation method for each signal
    end
    
    methods
        
        function this = var_step_signal_gen(signals, cp, method, p0)
            
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
                imethod = method;
                method  = cell(1,numel(signals));
                for is = 1:numel(signals)
                    method{is} = imethod;
                end
            end
            
            if ~exist('p0', 'var')
                p0 = [];
            end
            
            this.signals = signals;
            this.cp = cp;
            this.params = {};
            this.method = method;
            
            % fill in time steps params
            for kstep = 1:cp-1
                this.params= [this.params { ['dt_u' num2str(kstep-1)] }];
                this.p0(1,kstep) = 1;
            end
            
            for ku = 1:numel(signals)
                for k = 1:cp
                    this.params= [this.params { [signals{ku} '_u' num2str(k-1)]} ];
                end
            end
            if ~isempty(p0)
                this.p0 = p0;
            else
                this.p0 = zeros(numel(this.params),1);
            end
                       
           this.params_domain = repmat(BreachDomain(), 1, numel(this.params));
           this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));
            
        end
        
        function [X, time] = computeSignals(this,p, time) % compute the signals
            
            if size(p,1) ==1
                p = p';
            end
            
            cp = this.cp;
            X = zeros(numel(this.signals),numel(time));
            
            dt_cp = p(1:cp-1);
            t_cp = unique( [0; cumsum(dt_cp)]);
            
            pts_x = p(cp:end);
            for i_sig = 1:numel(this.signals)
                x_values = pts_x(1:cp);
                meth = this.method{i_sig};
                pts_x = pts_x(cp+1:end);
                
                if numel(t_cp)==1
                    x = x_values(end)*ones(numel(time),1);
                else
                    x_values = x_values(1:min([numel(t_cp) numel(x_values)]));
                    t_cp = t_cp(1:min([numel(t_cp) numel(x_values)]));
                    x = interp1(t_cp, x_values, time', meth, 'extrap');
                end
                X(i_sig,:) = x';
            end
        end
        
        function type = getType(this)
            type = 'varstep';
        end
    end
end
