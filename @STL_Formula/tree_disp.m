function [st, st_ids, not_expanded] =  tree_disp(phi,l, l_max)
    phis= get_children(phi);
    prefix= ''; 
    if l>0
        prefix = [repmat('  ', 1, l-1) '|-'];        
    end
    if nargin<3
        l_max =inf;
    end
    
    st_id = get_id(phi);
    if isempty(phis) 
        st = {[prefix disp(phi,0)]};                
        st_ids = {st_id};
        not_expanded = {[]}; 
    elseif numel(phis)==1
        [st1, st_ids1, not_expanded1] = tree_disp(phis{1}, l+1, l_max);        
        st = [ {[prefix disp(phi,2)]}; st1];        
        st_ids = [ st_id; st_ids1];        
        not_expanded = [{[]};not_expanded1];
    elseif l>l_max
        st = {[prefix st_id]} ;
        st_ids = {st_id};
        not_expanded={phi};
    else
        [st1, st_ids1, not_expanded1] = tree_disp(phis{1}, l+1, l_max);        
        [st2, st_ids2, not_expanded2] = tree_disp(phis{2}, l+1, l_max);        
        
        st = [ {[prefix disp(phi,2)]}; st1;st2];
        st_ids = [ st_id; st_ids1;st_ids2];
        not_expanded = [{[]} ;not_expanded1;not_expanded2];
    end
    
end