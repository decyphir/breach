function phi = set_in_signal_names(phi, in_signal_names)
%SET_IN_SIGNAL_NAMES defines input signals in a formula
%
%  syntax : phi = set_in_signal_names(phi, in_signal_names )
%
%  phi = set_in_signal_names(phi, in_signal_names) assumes in_signal_names
%  is an array of strings of the form sig_name_1, ..., sig_name_n 
%  where each sig_name_i is a name of an input signal

global BreachGlobOpt
    used_signals = STL_ExtractSignals (phi);
    
    % Add only input signals that do appear in the formula
    for(i=1:length(in_signal_names))
        for(j=1:length(used_signals))
            name1 = in_signal_names{i};
            name2 = used_signals{j};
            flag = strcmp(name1, name2);
            if (flag == 1)
                phi.in_signal_names{end+1} = name1;
                break;
            end
        end
    end
    
    % Propagate recursively the input signals to sub-formulas
    switch(phi.type)
    
        case 'predicate'
        case 'not'
            phi.phi = set_in_signal_names(phi.phi, phi.in_signal_names);
        case 'or'
            phi.phi1 = set_in_signal_names(phi.phi1, phi.in_signal_names);
            phi.phi2 = set_in_signal_names(phi.phi2, phi.in_signal_names);        
        case 'and'
            phi.phi1 = set_in_signal_names(phi.phi1, phi.in_signal_names);
            phi.phi2 = set_in_signal_names(phi.phi2, phi.in_signal_names);
        case 'andn'
            n_phi = numel(phi.phin);
            for i=1:n_phi
                phi.phin(i) = set_in_signal_names(phi.phin(i), phi.in_signal_names);
            end
        case '=>'
            phi.phi1 = set_in_signal_names(phi.phi1, phi.in_signal_names);
            phi.phi2 = set_in_signal_names(phi.phi2, phi.in_signal_names);
        case 'always'
            phi.phi = set_in_signal_names(phi.phi, phi.in_signal_names);
        case 'av_eventually'
            phi.phi = set_in_signal_names(phi.phi, phi.in_signal_names);
        case 'eventually'
            phi.phi = set_in_signal_names(phi.phi, phi.in_signal_names);
        case 'until'
            phi.phi1 = set_in_signal_names(phi.phi1, phi.in_signal_names);
            phi.phi2 = set_in_signal_names(phi.phi2, phi.in_signal_names);
    end
    
    % make sure the base formula gets updated with new parameters
    BreachGlobOpt.STLDB(phi.id) = phi;
end