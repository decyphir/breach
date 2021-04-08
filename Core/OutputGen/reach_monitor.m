classdef reach_monitor < req_monitor
% reach_monitor works for one signal at a time. Maintains a timed enveloppe
% and return the distance to it when outside to use for Falsification.
    
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
            
            I = [0 time(2)]; % Assumes constant delta for this time 
            
            [t_top, val_top] = RobustEv(t, Xin, I);
            Xin_top1 = interp1(t_top, val_top, time, 'previous', 'extrap');                                    
            Xin_top2 = [Xin_top1(1) Xin_top1(1:end-1)];
            Xin_top = max([Xin_top1;Xin_top2]);
            top = max([Xin_top;old_top]);
            
            [t_bot, val_bot] = RobustEv(t, -Xin, I);
            val_bot = -val_bot;
            Xin_bot1 = interp1(t_bot, val_bot, time, 'previous', 'extrap');                                    
            Xin_bot2 = [Xin_bot1(1) Xin_bot1(1:end-1)];                         
            Xin_bot = max([Xin_bot1;Xin_bot2]);
            
            bot = min([Xin_bot;old_bot]);
            
            m(this.signals_in{1}) = [time; top; bot];
            
            penalty_inside = min(min(old_top - Xin_top) , min(Xin_bot - old_bot))+... % this should force optimizer to get closer to one boundary
                             1e-9/(1+(sum(old_top==Xin_top)+sum(Xin_bot == old_bot)));    % Avoid zero if boundary is reached 
            
            if penalty_inside>=0
                val = penalty_inside;
            else
                val = -sum(top-bot); % not sure this is actually useful (-1 
                                     % could just work the same for the purpose 
                                     % of signaling we found a new bound,
                                     % though this provides info on the
                                     % volume of the new enveloppe. Somehow.               
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
       
        function f= plot_enveloppe(this, ax, time)
            m = this.RangesMap;
            if nargin<3
                ranges = m(this.signals_in{1});
                time= ranges(1,:);
            end
            if nargin<2
                ax = gca;                
            end
            
            axes(ax);
            hold on;
            [~, env] = this.computeSignals(time,0,0);
            vbot = env(2,:);
            vtop = env(1,:);
            
            plot(time, vbot, 'k', 'LineWidth', .5);
            plot(time, vtop, 'k', 'LineWidth', .5);
            t2 = [time, fliplr(time)];
            inBetween = [vbot, fliplr(vtop)];
            f = fill(t2, inBetween, 'k', 'FaceAlpha', 0.1);                        
        end        
        
    end
  
        
end
       
