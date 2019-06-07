function phi = set_out_signal_names(phi, out_signal_names)
%SET_OUT_SIGNAL_NAMES defines output signals in a formula
%
%  syntax : phi = set_out_signal_names(phi, out_signal_names )
%
%  phi = set_out_signal_names(phi, out_signal_names) assumes 
%  out_signal_names is an array of strings of the form 
%  sig_name_1, ..., sig_name_n 
%  where each sig_name_i is a name of an output signal.

global BreachGlobOpt
    used_signals = STL_ExtractSignals (phi);
    
    % Add only output signals that do appear in the formula
    for(i=1:length(out_signal_names))
        for(j=1:length(used_signals))
            name1 = out_signal_names{i};
            name2 = used_signals{j};
            flag = strcmp(name1, name2);
            if (flag == 1)
                phi.out_signal_names{end+1} = name1;
                    break;
            end
        end
    end
    
    % Propagate recursively the output signals to sub-formulas
    switch(phi.type)
    
        case 'predicate'
        case 'not'
            phi.phi = set_out_signal_names(phi.phi, phi.out_signal_names);
        case 'or'
            phi.phi1 = set_out_signal_names(phi.phi1, phi.out_signal_names);
            phi.phi2 = set_out_signal_names(phi.phi2, phi.out_signal_names);        
        case 'and'
            phi.phi1 = set_out_signal_names(phi.phi1, phi.out_signal_names);
            phi.phi2 = set_out_signal_names(phi.phi2, phi.out_signal_names);
        case 'andn'
            n_phi = numel(phi.phin);
            for i=1:n_phi
                phi.phin(i) = set_out_signal_names(phi.phin(i), phi.out_signal_names);
            end
        case '=>'
            phi.phi1 = set_out_signal_names(phi.phi1, phi.out_signal_names);
            phi.phi2 = set_out_signal_names(phi.phi2, phi.out_signal_names);
        case 'always'
            phi.phi = set_out_signal_names(phi.phi, phi.out_signal_names);
        case 'av_eventually'
            phi.phi = set_out_signal_names(phi.phi, phi.out_signal_names);
        case 'eventually'
            phi.phi = set_out_signal_names(phi.phi, phi.out_signal_names);
        case 'until'
            phi.phi1 = set_out_signal_names(phi.phi1, phi.out_signal_names);
            phi.phi2 = set_out_signal_names(phi.phi2, phi.out_signal_names);
    end
    
    % make sure the base formula gets updated with new parameters
    BreachGlobOpt.STLDB(phi.id) = phi;
end