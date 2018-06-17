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
            ax1 = F.AddAxes();
            F.AddSignals(signals_pre, ax1);
             
            ax2 = F.AddAxes();
            F.AddSignals(setdiff(this.signals_in, signals_pre), ax2);
            sig= this.signals{end};
            F.HighlightFalse(sig,ax2);
        end
    
        
    end
    
end