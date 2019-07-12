classdef alw_monitor < stl_monitor
    properties
        interval
        subphi
    end
    
    methods
        function this = alw_monitor(formula)
            this = this@stl_monitor(formula);
            this.subphi = get_children(this.formula);
            this.subphi = this.subphi{1};
            this.signals = {[this.formula_id '_violation'],... 
                            [this.formula_id '_quant_sat']};
                
            this.interval = get_interval(formula);
            this.init_P();
        end
        
        function [time, Xout] = computeSignals(this, time, X, p)
            this.init_tXp(time,X,p);

            % compute robustnes of top formula
            [time, Xout] = this.get_standard_rob(this.subphi, time);
            idx  = this.get_time_idx_interval(time,p);
            Xout =  [Xout(1,:)<0 ;...  % violation flags
                     Xout];         
            Xout(end-1:end, ~idx) = NaN;
   
        end

        
        function [v, t, Xout] = eval(this, t, X,p)
            [t, Xout] = this.computeSignals(t, X,p);
            v = min(Xout(end,:));
        end
        
        function st = disp(this)
            phi = this.subphi;
            st = sprintf(['%s := alw_%s (%s)\n'], this.formula_id,this.interval,  disp(phi,1));
            
            if ~strcmp(get_type(phi),'predicate')
                st = [st '  where\n'];
                predicates = STL_ExtractPredicates(phi);
                for ip = 1:numel(predicates)
                    st =   sprintf([ st '%s := %s \n' ], get_id(predicates(ip)), disp(predicates(ip)));
                end
            end
            
            if nargout == 0
                fprintf(st);
            end
     
        end
        
        function plot_diagnosis(this, F)
            % Assumes F has data about this formula 
            
            F.AddSignals(this.signals_in);
            sig= this.signals{end};
            F.HighlightFalse(sig);
        
        end
    
    end
    
    methods  (Access=protected)
        function idx__ = get_time_idx_interval(this, t__, p__)
            this.assign_params(p__);
            interval = eval(this.interval);
            idx__ = (t__>= interval(1))&(t__<=interval(end));
            if ~any(idx__)
                idx__(end) = 1;
            end
        end
        
        
    end
end