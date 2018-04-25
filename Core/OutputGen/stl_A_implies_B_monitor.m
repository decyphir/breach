classdef stl_A_implies_B_monitor < stl_monitor
    properties
    pre
    post
    end
    methods
        function this = stl_A_implies_B_monitor(formula)
            this = this@stl_monitor(formula);
            
            if ~isequal(get_type(this.formula), '=>')
                error('formula is not of type A => B');
            end

            subphis = get_children(formula);
            this.pre  = subphis{1};
            this.post = subphis{2}; 
            
        end
        
        function h= plot_diagnosis(this, time, X, p )
            [tau, val] = this.computeSignals(time, X, p);
            [valA, tauA] = STL_Eval(this.Sys, this.pre, this.P,this.P.traj{1}, tau);
            [valB, tauB] = STL_Eval(this.Sys, this.post, this.P,this.P.traj{1}, tau);
       
            ax(1) = subplot(3,1,1);
            stairs(tauA, valA);
            grid on;
            title(sprintf('pre: %s', disp(this.pre)))
           highlight_truth_intervals(tauA,valA);
            
            ax(2) = subplot(3,1,2);
            stairs(tauB, valB)
            grid on;
            title(sprintf('post: %s', disp(this.post)))
            highlight_truth_intervals(tauB,valB);
            
            ax(3) =subplot(3,1,3);
            stairs(tau, val);
            grid on;
            title(sprintf('%s', disp(this.formula)))
           highlight_truth_intervals(tau,val);
           highlight_truth_intervals(tauA,valA, 'g', 0, [0.6 0 0.6], 0.3);
            linkaxes(ax,'x');
            
        end
    end
    
    
end