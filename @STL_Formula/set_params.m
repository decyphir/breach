function phi = set_params(phi, P)
%SET_PARAMS adds a parameter structure to a formula
%  
if ~isempty(fieldnames(P))
    phi.params.default_params = P; 
    % make sure the base formula gets updated with new parameters     
    assignin('base', phi.id, phi);
end
