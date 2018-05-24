function subphis = get_children(phi)

subphis = [];
if ~isempty(phi.phi) 
    subphis = {phi.phi};
elseif ~isempty(phi.phi1)
    subphis = {phi.phi1, phi.phi2};
end