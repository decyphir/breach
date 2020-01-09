classdef reach_monitor < req_monitor

    
    properties (Constant)
        RangesMap = containers.Map()
    end
    
    
    methods
        function this = reach_monitor(name, signal, time)
            this.name = name;
            top_name = [signal '_top_'];
            bot_name = [ signal '_bot_'];                                   
            
            this.signals = {top_name bot_name};                        
            this.signals_in = {signal};
            
            m = this.RangesMap;
            m(signal) = [time;  0*time-inf; 0*time+inf];                        
        
        end
        
            
        function [t, Xout] = computeSignals(this, t, Xin, pin)
            % computes distance to bounds top and bot 
            
            ranges = this.RangesMap(this.signals_in{1});
            time= ranges(1,:);
            top = ranges(2,:);
            bot = ranges(3,:);
            
            Xout = [interp1(time, top, t, 'linear', 'extrap'); ...
                    interp1(time, bot, t, 'linear', 'extrap')];                        
        
        end
    
        function [val, t, Xout] = eval(this,t, Xin, pin)
           
            m = this.RangesMap;            
            ranges = m(this.signals_in{1});
            time= ranges(1,:);
            old_top = ranges(2,:);
            old_bot = ranges(3,:);
            
                        
            % Update top and bot
            Xin_t = interp1(t, Xin, time, 'linear', 'extrap');                                   
                        
            top = max([Xin_t ;old_top]);
            bot = min([Xin_t ;old_bot]);            
            
            m(this.signals_in{1}) = [time; top; bot];
            
            penalty_inside = min(min(old_top - Xin_t) , min(Xin_t - old_bot));
            
            if penalty_inside>0
                val = penalty_inside;
            else
                val = -sum(top-bot);
            end
            
            [t, Xout] = this.computeSignals(t, Xin, pin);
        end
           
 

        function varargout = disp(this)
            st= [this.name ': Reach monitor for signal ' this.signals_in{1} '.\n'];
            
            
            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end
            
        end
        
    end
  
        
end
       
