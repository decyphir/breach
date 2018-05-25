function params = get_params_from_attributes(S, att)
fname = [att 's_idx'];
if ~isfield(S, fname)
    error('%s not a param attribute', att);
end
idx = S.(fname);

params = S.params(idx);
