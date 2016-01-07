classdef fixed_cp_signal_gen < signal_gen
    
    % Control points, uni/fixed step signal generation 
    properties 
        cp       % number of control points for each signal
        method   % interpolation method for each signal
    end
    
    methods 
        
        function this = fixed_cp_signal_gen(signals, cp, method)
           this.signals = signals;
           this.cp = cp;
           this.params = {};
           this.method = method;
           for ku = 1:numel(signals)
               for k = 1:cp(ku)
                   this.params= {this.params{:} [signals{ku} '_u' num2str(k-1)]};
               end
           end
           this.p0= zeros(numel(signals)*sum(cp),1);           
        end
            
        function X = computeSignals(this,p, time) % compute the signals
            if size(p,1) ==1
                p = p';
            end

            X = zeros(numel(this.cp),numel(time));
            pts_x = p; 
            for i_cp = 1:numel(this.cp)
                cp_values = pts_x(1:this.cp(i_cp));
                meth = this.method{i_cp};
                pts_x = pts_x(this.cp(i_cp)+1:end);
                t_cp = linspace(time(1), time(end), this.cp(i_cp)+1)';
                if numel(t_cp)==2
                    x = cp_values(1)*ones(numel(time),1);
                else
                    x = interp1(t_cp(1:end-1), cp_values, time', meth, 'extrap');
                end
                X(i_cp,:) = x';
            end        
        end
        
        function type = getType(this)
            type = 'unistep';
        end
    end    
end
