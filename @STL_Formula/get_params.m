function params = get_params(phi)

if isfield(phi.params,'default_params')
    params = phi.params.default_params;
else
    params = struct;
end
end
