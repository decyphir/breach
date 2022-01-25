classdef var_cp_signal_gen < signal_gen
    % var_cp_signal_gen  A class derived from signal_gen to generate signals from control points
    %                    with variable time steps, using different interpolation methods.
    %
    % var_cp_signal_gen Properties
    %   cp     - number of control points
    %   method - interpolation method for each signal, accepts all methods provided by interp1
    %
    % var_cp_signal_gen Methods
    %   var_cp_signal_gen -  constructor, takes signal names, number of
    %                        control points, interpolation methods and a
    %                        p0.
    %  See also signal_gen.
    
    
    properties
        num_cp       % number of control points for each signal
        method   % interpolation method for each signal
    end
    
    methods
        
        function this = var_cp_signal_gen(signals, cp, method, p0)
            
            if nargin == 1
                cp = 1;
                method = 'previous';
            elseif nargin == 2
                method = 'previous';
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
            this.num_cp = cp;
            this.params = {};
            this.method = method;
            this.p0 = p0;
            
            p0 = zeros(2*sum(cp)-numel(signals),1);
            for ku = 1:numel(signals)
                k=0;                
                for k = 1:cp(ku)-1
                    this.params= [this.params { [signals{ku} '_u' num2str(k-1)] [signals{ku} '_dt' num2str(k-1)] }];
                    if isempty(this.p0)
                        p0(numel(this.params))=1;
                    end
                end
                if isempty(k) % case where there is only one control point 
                  k=0;
                end
                this.params = [this.params [signals{ku} '_u' num2str(k)]];
            end
            
            if isempty(this.p0)
                this.p0=p0;
            end
        
           this.params_domain = repmat(BreachDomain(), 1, numel(this.params));
           this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));
                 
        end
                        
        function [cp_times, cp_values, t_bnd, cp_bnd] = get_cp(this, signal)
            
            i_sig = find(strcmp(signal,this.signals),1);
            i_p = sum(this.num_cp(1:i_sig-1))+1;
            idx_cp = i_p:2:i_p+2*this.num_cp(i_sig)-2;
            idx_t = i_p+1:2:i_p+2*this.num_cp(i_sig)-2;
            
            cp_values = this.p0(idx_cp);
            cp_dt  = this.p0(idx_t); 
            
            cp_times = cp_values*0;
            for ic = 1:numel(cp_dt)
               cp_times(ic+1) = sum(cp_dt(1:ic));
            end            
                                   
            cp_bnd = [cp_values cp_values];
            t_bnd = [cp_times cp_times];
                        
            for ic = 1:numel(idx_cp)
                dom_cp(ic) = this.params_domain(idx_cp(ic));
                if ~isempty(dom_cp(ic).domain)
                    cp_bnd(ic,1) = dom_cp(ic).domain(1);
                    cp_bnd(ic,2) = dom_cp(ic).domain(2);                    
                end
                if ic>1
                   dom = this.params_domain(idx_t(ic-1));
                   if ~isempty(dom.domain)                                              
                       t_bnd(ic,1) = sum(t_bnd(1:ic-1,1))+dom.domain(1);
                       t_bnd(ic,2) = sum(t_bnd(1:ic-1,2))+dom.domain(2);
                   else
                       t_bnd(ic,1) = sum(t_bnd(1:ic,1));
                       t_bnd(ic,2) = sum(t_bnd(1:ic,2));                       
                   end
                end                
            end
        end
        
                                        
        function [X, time] = computeSignals(this,p, time) % compute the signals
            
            if size(p,1) ==1
                p = p';
            end
            
            
            X = zeros(numel(this.num_cp),numel(time));
            pts_x = p;
            for i_cp = 1:numel(this.num_cp)
                cp_values = pts_x(1:2*this.num_cp(i_cp)-1);
                meth = this.method{i_cp};
                pts_x = pts_x(2*this.num_cp(i_cp):end);
                
                dt_cp = cp_values(2:2:end-1);
                x_values = cp_values(1:2:end);
                dt_cp(dt_cp==0) = 1e-10;
                t_cp = [0; cumsum(dt_cp)];
                [t_cp, i_t_cp] = unique( [0; cumsum(dt_cp)], 'last');                
                x_values = x_values(i_t_cp);

                if numel(t_cp)==1
                    x = x_values(end)*ones(numel(time),1);
                else
                    x_values = x_values(1:min([numel(t_cp) numel(x_values)]));
                    t_cp = t_cp(1:min([numel(t_cp) numel(x_values)]));
                    x = interp1(t_cp, x_values, time', meth, 'extrap');
                end
                X(i_cp,:) = x';
            end
        end
                
        function plot(this, signal, time)
            
            [t_cp, x_cp] = this.get_cp(signal);
            time =  sort([time t_cp']);
            plot@signal_gen(this,signal, time);
            
            % plot control points
            hold on;
            plot(t_cp,x_cp,'or','MarkerSize', 6, 'MarkerFaceColor', [1 0 0]);
            legend({signal,'control points (cp)'}, 'Interpreter', 'None');
            
        end                        
          
        function type = getType(this)
            type = 'varstep';
        end
        
        function args = getSignalGenArgs(this)
            args = {'num_cp','method'};         
        end
        
    end
end
