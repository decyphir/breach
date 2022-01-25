classdef cp_signal_gen < signal_gen
    % cp_signal_gen  A class derived from signal_gen to generate signals from control points using different interpolation methods.
    %
    % cp_signal_gen Properties
    %   cp     - number of control points
    %   method - interpolation method for each signal, accepts all methods provided by interp1
    %
    % disturb_cp_signal_gen Methods
    %   cp_signal_gen -  constructor, takes signal names, number of
    %                            control points, interpolation methods and
    %                            ref control points p0.
    %  See also signal_gen.
    
    
    properties
        cp       % number of control points for each signal
        method   % interpolation method for each signal
    end
    
    methods
        
        function this = cp_signal_gen(signals, cp, method, p0)
            
            if nargin == 1
                cp = 1;
                method = 'linear';
            elseif nargin == 2
                method = 'linear';
            end
            
            if ~iscell(signals)
                signals = {signals};
            end
            
            if ~iscell(method)
                imethod = method;
                method  = cell(1,numel(signals));
                for is = 1:numel(signals)
                    method{is} = imethod;
                end
            end
            
            if ~exist('p0', 'var')
                p0=[];
            end
            
            this.signals = signals;
            this.cp = cp;
            this.params = {};
            this.method = method;
            this.p0 = p0;
            
            p0 = zeros(2*sum(cp)-numel(signals),1);
            for ku = 1:numel(signals)
                this.params = [this.params [ signals{ku} '_u0']];
                for k = 1:cp(ku)-1
                    this.params= [this.params {  [signals{ku} '_t' num2str(k)] [signals{ku} '_u' num2str(k)]}];
                    if isempty(this.p0)
                        p0(numel(this.params)-1)=k;
                    end
                end
            end
            
            if isempty(this.p0)
                this.p0=p0;
            end
            
            
           this.params_domain = repmat(BreachDomain(), 1, numel(this.params));
           this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));
 
        end
        
        function [X, time] = computeSignals(this,p, time) % compute the signals
            
            if size(p,1) ==1
                p = p';
            end
            
            X = zeros(numel(this.cp),numel(time));
            pts_x = p;
            for i_cp = 1:numel(this.cp)
                cp_values = pts_x(1:2*this.cp(i_cp)-1);
                meth = this.method{i_cp};
                pts_x = pts_x(2*this.cp(i_cp):end);
                
                t_cp = [0; cp_values(2:2:end-1)];
                x_values = cp_values(1:2:end);
                if numel(t_cp)==1
                    x = x_values(end)*ones(numel(time),1);
                else
                    x_values = x_values(1:min([numel(t_cp) numel(x_values)]));
                    t_cp = t_cp(1:min([numel(t_cp) numel(x_values)]));
                    % TODO sort t_cp
                    [t_cp, order] = sort(t_cp);
                    x_values = x_values(order);
                    x = interp1(t_cp, x_values, time', meth, 'extrap');
                end
                X(i_cp,:) = x';
            end
        end
        
        function type = getType(this)
            type = 'varstep';
        end
    end
end
