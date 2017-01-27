function [phi, P] = check_params(phi)

if isfield(phi.params,'default_params')
    P = phi.params.default_params;
else
    P = struct;
end

if ~isempty(phi.phi)
    [~, P0] = check_params(phi.phi);
    fnames =  fieldnames(P0);
    for ifn = 1:numel(fnames)
        P.(fnames{ifn}) = P0.(fnames{ifn});
    end
end

if ~isempty(phi.phi1)
    [~, P0] = check_params(phi.phi1);
    fnames =  fieldnames(P0);
    for ifn = 1:numel(fnames)
        P.(fnames{ifn}) = P0.(fnames{ifn});
    end
end

if ~isempty(phi.phi2)
    [~, P0] = check_params(phi.phi2);
    fnames =  fieldnames(P0);
    for ifn = 1:numel(fnames)
        P.(fnames{ifn}) = P0.(fnames{ifn});
    end
end

phi.params.default_params=P;

end
