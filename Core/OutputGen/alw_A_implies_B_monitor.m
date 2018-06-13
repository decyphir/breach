classdef alw_A_implies_B_monitor < alw_monitor
    properties
        pre
        post
    end
    methods
        function this = alw_A_implies_B_monitor(formula)
            this = this@alw_monitor(formula);
          
            subphis = get_children(formula);
            if ~isequal(get_type(subphis{1}), '=>')
                error('formula is not of type alw A => B');
            end
            
            this.pre  = subphis{1};
            this.post = subphis{2}; 
            
        end
        
        function plot_diagnosis(this, F)
            % Assumes F has data about this formula 
            
            F.AddSignals(this.signals_in);
            sig= this.signals{end};
            F.HighlightFalse(sig);
        end
    
        
    end
    
end