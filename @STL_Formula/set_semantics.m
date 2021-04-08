function phi = set_semantics(phi, semantics, set_all)
% SET_SEMANTICS assign a new robust semantics to a formula

% If set_all is true, we will set the semantics of all underlying
% subformulas as well
% If set_all is false, ONLY this formula (no subformulas) will have the new
% semantics.

if nargin < 3
    % Standard case: set_all is true 
    set_all = 1;
end

if set_all
    % Set new semantics for ALL subformulas
    switch (phi.type)
        case 'predicate'
            % No subformula to set semantics for
            
        case 'not'
            phi.phi = set_semantics(phi.phi, semantics);
            
        case{'always', 'eventually', 'historically', 'once'}
            phi.phi = set_semantics(phi.phi, semantics);
            
        case {'and', 'or', '=>'}
            phi.phi1 = set_semantics(phi.phi1, semantics);
            phi.phi2 = set_semantics(phi.phi2, semantics);
            
        case 'until'
            phi.phi1 = set_semantics(phi.phi1, semantics);
            phi.phi2 = set_semantics(phi.phi2, semantics);
    end
end

% Only set new semantics for the current formula
phi.semantics = semantics;

global BreachGlobOpt
BreachGlobOpt.STLDB(phi.id) = phi;