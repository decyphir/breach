function phi = set_in_signal_names(phi, in_signal_names,reset)
%SET_IN_SIGNAL_NAMES defines input signals in a formula
%
%  syntax : phi = set_in_signal_names(phi, in_signal_names )
%
%  phi = set_in_signal_names(phi, in_signal_names) assumes in_signal_names
%  is an array of strings of the form sig_name_1, ..., sig_name_n 
%  where each sig_name_i is a name of an input signal

global BreachGlobOpt

if nargin<3
    reset = 0;
end

switch reset 
    case 0  % checks and add
        % Add only output signals that do appear in the formula
        used_signals = STL_ExtractSignals (phi);
        new_in_signals = intersect(used_signals, in_signal_names);        
        phi.in_signal_names = union(phi.in_signal_names, new_in_signals);
    
    case 1  % checks and reset         
        % Add only output signals that do appear in the formula
        used_signals = STL_ExtractSignals (phi);        
        phi.in_signal_names = intersect(used_signals, in_signal_names);        
                 
    case 2 % don't check and reset
        phi.in_signal_names = in_signal_names;

end

% Propagate recursively the output signals to sub-formulas
switch(phi.type)
    
    case 'predicate'
    case 'not'
        phi.phi = set_in_signal_names(phi.phi, phi.in_signal_names, 2);
    case 'or'
        phi.phi1 = set_in_signal_names(phi.phi1, phi.in_signal_names, 2);
        phi.phi2 = set_in_signal_names(phi.phi2, phi.in_signal_names, 2);
    case 'and'
        phi.phi1 = set_in_signal_names(phi.phi1, phi.in_signal_names, 2);
        phi.phi2 = set_in_signal_names(phi.phi2, phi.in_signal_names, 2);
    case 'andn'
        n_phi = numel(phi.phin);
        for i=1:n_phi
            phi.phin(i) = set_in_signal_names(phi.phin(i), phi.in_signal_names, 2);
        end
    case '=>'
        phi.phi1 = set_in_signal_names(phi.phi1, phi.in_signal_names, 2);
        phi.phi2 = set_in_signal_names(phi.phi2, phi.in_signal_names, 2);
    case 'always'
        phi.phi = set_in_signal_names(phi.phi, phi.in_signal_names, 2);
    case 'av_eventually'
        phi.phi = set_in_signal_names(phi.phi, phi.in_signal_names, 2);
    case 'eventually'
        phi.phi = set_in_signal_names(phi.phi, phi.in_signal_names, 2);
    case 'until'
        phi.phi1 = set_in_signal_names(phi.phi1, phi.in_signal_names, 2);
        phi.phi2 = set_in_signal_names(phi.phi2, phi.in_signal_names, 2);
end

    
    % make sure the base formula gets updated with new parameters
    BreachGlobOpt.STLDB(phi.id) = phi;
end