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
        
        
        
        function ax= plot_diagnostics(this, time, X, p, style)
            [tau, val] = this.computeSignals(time, X, p);
            [valA, tauA] = STL_Eval(this.Sys, this.pre, this.P,this.P.traj{1}, tau);
            [valB, tauB] = STL_Eval(this.Sys, this.post, this.P,this.P.traj{1}, tau);
            
            if ~exist('style', 'var')
                style = '';
            end
            
            if ~strcmp(style, 'compact')
                ax(1) = subplot(3,1,1);
                stairs(tauA, valA);
                grid on;
                title(sprintf('pre: %s', disp(this.pre)), 'Interpreter', 'None');
                highlight_truth_intervals(tauA,valA);
                
                ax(2) = subplot(3,1,2);
                stairs(tauB, valB)
                grid on;
                title(sprintf('post: %s', disp(this.post)), 'Interpreter', 'None');
                highlight_truth_intervals(tauB,valB);
                ax(3) =subplot(3,1,3);
                linkaxes(ax,'x');
            else
                %stairs(tau, val);
                grid on;
                title(sprintf('%s', get_id(this.formula)), 'Interpreter', 'None');
                highlight_truth_intervals(tau,val);
                highlight_truth_intervals(tauA,valA, 'g', 0, [0.6 0 0.6], 0.3);
                ax = gca;
            end
            
        end
    end
    
    
end