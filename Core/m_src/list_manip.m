classdef list_manip
    methods (Static)
        
        function l = append_last(l, elem)
            if ~iscell(elem)
                elem = {elem};
            end
            l = [l elem];
        end
        
        function l = append_first(l, elem)
            if ~iscell(elem)
                elem = {elem};
            end
            l = [elem l];
        end
        
        function l = insert_elem(l,elem,i)
            n = numel(l);
            if n<i
                l{i}=elem;
            else
                l = [l(1:i-1) {elem} l(i:end)];
            end
        end
        
        function idx = find_elem(l,elem)
            idx = find(cellfun(@(c)(isequal(c,elem)),l),1);
        end
        
        function idx = find_all_elem(l,elem)
            idx = find(cellfun(@(c)(isequal(c,elem)),l));
        end
                
        function idx = find_str(l,elem)
            idx = find(strcmp(l, elem),1);
        end
        
        function l = move_up(l,i)
            n = numel(l);
            if i==2
                l = [l(i) l(1) l(i+1:end)];
            elseif i>2&&i<=n
                l = [l(1:i-2) l(i) l(i-1) l(i+1:end)];
            end            
        end
        
        function l = move_down(l,i)
            n = numel(l);
            if i==2
                l = [l(i) l(1) l(i+1:end)];
            elseif i>2&&i<=n
                l = [l(1:i-2) l(i) l(i-1) l(i+1:end)];
            end            
        end
        
        function l = remove_elem(l,i)
            n=numel(l);
            if i==n
                l = l(1:i-1);
            elseif i==1
                l = l(2:end);
            elseif i>1&&i<n
                l = [l(1:i-1) l(i+1:end)];
            end            
        end
        
        function l = remove_empty(l)    
            idx_not_empty = cellfun(@(c)(~isempty(c)), l); 
            l = l(idx_not_empty);            
        end
        
        function st =  to_string(l, sep)
            st = l{1};
            if nargin <2
                sep= ',';
            end    
            if numel(l)>1
                for n = 2:numel(l)
                    st=  [st sep l{n}];
                end
            end            
        end
    end
    
end