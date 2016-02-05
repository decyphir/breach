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
        cp       % number of control points for each signal
        method   % interpolation method for each signal
    end
    
    methods 
        
        function this = var_cp_signal_gen(signals, cp, method)
           this.signals = signals;
           this.cp = cp;
           this.params = {};
           this.method = method;
           for ku = 1:numel(signals)
               for k = 1:cp(ku)
                   this.params= {this.params{:} [signals{ku} '_dt' num2str(k-1)] [signals{ku} '_u' num2str(k-1)]};
               end
           end
           this.p0= zeros(2*sum(cp)*numel(signals),1);
           this.p0(1:2:end) = 1;
        end
        
        function X = computeSignals(this,p, time) % compute the signals
           
            if size(p,1) ==1
                p = p';
            end

            
            X = zeros(numel(this.cp),numel(time));
            pts_x = p;
            for i_cp = 1:numel(this.cp)
                cp_values = pts_x(1:2*this.cp(i_cp));
                meth = this.method{i_cp};
                pts_x = pts_x(2*this.cp(i_cp)+1:end);
                
                dt_cp = cp_values(1:2:end);
                t_cp = unique( [0; cumsum(dt_cp)]);
                x_values = cp_values(2:2:end);
                x_values = x_values(1:min([numel(t_cp) numel(x_values)]));
                t_cp = t_cp(1:min([numel(t_cp) numel(x_values)]));
                if numel(t_cp)==1
                    x = x_values(1)*ones(numel(time),1);
                else
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
