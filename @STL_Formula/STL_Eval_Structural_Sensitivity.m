function is_sensitive = STL_Eval_Structural_Sensitivity(Sys, phi, P, traj1, traj2, inout, relabs, taus)


%% default parameters
if isfield(phi.params,'default_params')
    pnames = fieldnames(phi.params.default_params);
    for ip = 1:numel(pnames)
        % if P doesn't define the parameter, use default
        if FindParam(P,pnames{ip})> size(P.pts,1)
            pval = phi.params.default_params.(pnames{ip});
            if isscalar(pval)
                P = SetParam(P,pnames{ip},pval );
            else
                P.(pnames{ip}) = pval;
            end
        end
    end
end

partition = [];
if (strcmp(inout, 'in'))
    partition = get_in_signal_names(phi);
elseif (strcmp(inout, 'out'))
    partition = get_out_signal_names(phi);
end

[trajs_sensitive, time_values__] = ...
    calculate_structural_sensitivity(Sys, phi, P, traj1, traj2, partition, relabs, taus);

% TODO: Can we ignore time values here?
is_sensitive = any(trajs_sensitive);


end
