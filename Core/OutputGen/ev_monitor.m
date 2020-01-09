classdef ev_monitor < alw_monitor
    
    methods
        
        function [v__, Xout__] = eval(this, t, X__,p__)
            [~, Xout__] = this.computeSignals(t, X__,p__);
            idx  = this.get_time_idx_interval(t,p__);
            v__ = max(Xout__(end,idx));
        end
        
        function varargout = disp(this)
            phi = this.formula;
            st = sprintf(['%s := ev_%s (%s)\n'], this.formula_id,this.interval,  disp(phi,1));
            
            if ~strcmp(get_type(phi),'predicate')
                st = [st '  where\n'];
                predicates = STL_ExtractPredicates(phi);
                for ip = 1:numel(predicates)
                    st =   sprintf([ st '%s := %s \n' ], get_id(predicates(ip)), disp(predicates(ip)));
                end
            end
            
            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end
            
        end
    end
 end