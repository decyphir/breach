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
               
        function [v, t, Xout] = eval(this, t, X,p)
            [~, Xout] = this.computeSignals(t, X,p);
            idx  = this.get_time_idx_interval(t,p);
            Xout(end-1,:) = Xout(end,:)<0;         % violation flags
            Xout(end-1:end, ~idx) = NaN;
            if isequal(get_type(this.post), 'eventually')
                int_post = this.get_post_interval(p);
                if int_post(2)<inf  % prob not necessary
                    a = int_post(1);
                    b = int_post(2);
                    s = Xout(end-1,:);
                    s(isnan(s))=0;
                    sev = 0*s;
                    Itrue = get_truth_intervals(s); % these are indexes 
                    for ii = 1:size(Itrue,1)
                        t1 = t(Itrue(ii,1));
                        t2 = t(Itrue(ii,2));
                        
                        tn1 = t1+a;
                        tn2 = t2+a+b;
                        
                        idx1 = max(Itrue(ii,1),find(t>=tn1,1)-1);
                        if tn2>t(end)
                            idx2 = numel(t);
                        else
                            idx2  = find(t>tn2,1)-1;
                        end
                        sev(idx1:idx2)=1; 
                    end
                    
                    Xout(end-1,:) = sev;
                end
            end
            v = min(Xout(end,idx));
        end
        
        function plot_diagnostics(this, F)
            % Assumes F has data about this formula 
            signals_pre = STL_ExtractSignals(this.pre);
            sig = this.signals{end-1};
            
            % plots pre signals
            ax1 = F.AddAxes();
            F.AddSignals(signals_pre, ax1);
            
            signals_post = setdiff(this.signals_in, signals_pre);
            if ~isempty(signals_post)
                % plots post intervals
                ax2 = F.AddAxes();
                F.AddSignals(signals_post, ax2);
                F.HighlightFalse(sig ,ax2, true);
            else
                F.HighlightFalse(sig ,ax1,true);
            end
        end
    
        function int_post__ = get_post_interval(this__,p__)
                this__.assign_params(p__);
                int__ = get_interval(this__.post);
                int_post__ =eval(int__);
        end
        
    end
    
end