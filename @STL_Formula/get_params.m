function params = get_params(phi)

if isfield(phi.params,'P')
    params = params.P;
else
    params = struct;
end
end
