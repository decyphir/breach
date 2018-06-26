classdef alw_A_implies_B_monitor < alw_monitor
    properties
        pre
        post
    end
    methods
        function this = alw_A_implies_B_monitor(formula)
            this = this@alw_monitor(formula);
          
            subphi = get_children(formula);
            if ~isequal(get_type(subphi{1}), '=>')
                error('formula is not of type alw A => B');
            end
            subphis = get_children(subphi{1});
            this.pre  = subphis{1};
            this.post = subphis{2}; 
            
        end
        
        function plot_diagnosis(this, F)
            % Assumes F has data about this formula 
            signals_pre = STL_ExtractSignals(this.pre);
            sig = this.signals{end};
            
            % plots pre intervals
            ax1 = F.AddAxes();
            F.AddSignals(signals_pre, ax1);
            int_false = F.HighlightFalse(sig,ax1);
            
            %  
            ax2 = F.AddAxes();
            F.AddSignals(setdiff(this.signals_in, signals_pre), ax2);
            
            % get intervals
            int_post = get_interval(this.post);
            if isequal(get_type(this.post), 'eventually')||~isempty(int_post) % the following works for A=>ev[t0, t1]B 
                itr= F.itraj;
                p = F.BrSet.GetParam(this.params, itr);
                this.assign_params(p);
                int_post =eval(int_post);
                if int_post(2)<inf
                    for ii = 1:size(int_false,1)
                        int_false_post(ii,: ) = [int_false(ii,1)+int_post(2) int_false(ii,2)+int_post(2)];
                        highlight_interval(ax2, int_false_post(ii,:), 'r', 0.3 );
                    end
                else
                    F.HighlightFalse(sig,ax2);
                end
            else
                F.HighlightFalse(sig,ax2);
            end
        end
    
        
    end
    
end