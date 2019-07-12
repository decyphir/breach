function map = get_formula_name_map(phi, map)
%GET_FORMULA_NAME_MAP creates a hash map of (formula_id, formula_text)
%pairs
%
%  syntax : map = get_formula_name_map(phi, map )
%
%  where phi is a STL_Formula and map is a Map of (string, string)

    id = get_id(phi);
    name = disp(phi,1);
    
    if(~map.isKey(id))
        map(id) = name;
    end
    
    switch(phi.type)
    
        case 'predicate'
            signals = STL_ExtractSignals(phi);
            for(i=1:length(signals))
                signal_name = signals{i};
                if(~map.isKey(signal_name))
                    map(signal_name) = signal_name;
                end
            end
        case 'not'
            map = get_formula_name_map(phi.phi, map);
        case 'or'
            map = get_formula_name_map(phi.phi1, map);
            map = get_formula_name_map(phi.phi2, map);
        case 'and'
            map = get_formula_name_map(phi.phi1, map);
            map = get_formula_name_map(phi.phi2, map);
        case 'andn'
            n_phi = numel(phi.phin);
            for i=1:n_phi
                map = get_formula_name_map(phi.phin(i), map);
            end
        case '=>'
            map = get_formula_name_map(phi.phi1, map);
            map = get_formula_name_map(phi.phi2, map);
        case 'always'
            map = get_formula_name_map(phi.phi, map);
        case 'av_eventually'
            map = get_formula_name_map(phi.phi, map);
        case 'eventually'
            map = get_formula_name_map(phi.phi, map);
        case 'until'
            map = get_formula_name_map(phi.phi1, map);
            map = get_formula_name_map(phi.phi2, map);
    end
end