function signals = get_signals_from_attributes(S, att)
fname = [att 's_idx'];
if ~isfield(S, fname)
    error('%s not a signal attribute', att);
end

idx = S.(fname);
signals = S.signals(idx);
