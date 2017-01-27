function [Pf, val_opt, nb_call_total] = Falsify (Sys, phi, falsif_opt, params_prop)
%FALSIFY Tries to falsify property phi
% 
% Synopsis: Pf = Falsify(Sys, phi, falsif_opt[, param_prop])
% 
% Input:
%  - Sys        : System Breach structure
%  - phi        : property to falsify
%  - falsif_opt : parameters for falsification - contains fields params, ranges, nb_init, nb_iter
%  - param_prop : parameters for phi -  contains fields params , values
% 
% Output:
%  - Pf : parameter set leading to falsification, or least robust
%         parameters found in case falsification failed
% 
% Example:
%  mdl = 'Autotrans_shift';
%  Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'UniStep1');
%  Sys.tspan = 0:.01:50;
%
%  falsif_prop.params_u = {'throttle_u0'};
%  falsif_opt.ranges = [0 100];
%
%  falsif_opt.nb_init = 10;
%  falsif_opt.nb_iter = 100;
%
%  Pr = Falsify(Sys, phi1, falsif_opt);
% 
%See also SOptimProp SOptimPropLog
%

%% Create system and input strategy

params = falsif_opt.params;
ranges = falsif_opt.ranges;
Pu = CreateParamSet(Sys, params, ranges);

if exist('params_prop','var')
    Pu = SetParam(Pu, params_prop.params,params_prop.values);
end

opt.timeout = getfield_default(falsif_opt, 'timeout', inf);
nb_restart = getfield_default(falsif_opt, 'nb_restart', 0); 
opt.nb_max_call = getfield_default(falsif_opt, 'nb_max_call', inf); 

nb_init = falsif_opt.nb_init;
nb_iter = falsif_opt.nb_iter;

%% Arguments for falsification

opt.OptimType = 'Min';
opt.MaxIter = nb_iter;
opt.StopWhenFoundInit = 1;
opt.StopWhenFound = 1;
opt.ubound = ranges(:,2);
opt.lbound = ranges(:,1);
opt.Ninit = 10; 
nb_init_total = 1;
nb_call_total = 0;

tic;
for restart = 0:nb_restart
    
    % Get a set of initial search points
    Pr_new = QuasiRefine(Pu, nb_init, nb_init_total);
    
    % Falsification step
    [val_opt, Propt, status,nb_call] = SOptimPropNM(Sys, Pr_new, phi, opt);
    ifalse = find(val_opt<0);
    tspent = toc;  

    if ~isempty(ifalse)
        Pf = Sselect(Propt, ifalse);
        fprintf(' ---- Falsified with\n');
        param_values = GetParam(Pf, params);
        for ip = 1:numel(params)
            fprintf( '        %s = %g\n', params{ip},param_values(ip))
        end
        fprintf('\n');
        return;
    else
        Pf = Propt;
        opt.timeout = opt.timeout - tspent;
        fprintf('Restart %d, time left: %g\n',restart, opt.timeout); 

        if opt.timeout <0
            fprintf('\nFalsification timed out.\n')
            fprintf(['Final robustness: ' num2str(val_opt)]);
            return;
        end
    
    end
    
    nb_init_total = nb_init_total + nb_init;
    nb_init = nb_init_total; % Fibo-increasing numbers of initials   
    nb_call_total = nb_call_total+nb_call;
    opt.nb_max_call = opt.nb_max_call - nb_call;
    
    % break if stop was requested
    if (status == -1)
        fprintf('\nFalsification stopped at user''s request.\n')       
        fprintf(['Final robustness: ' num2str(val_opt)]);
        return;    
    end
    % break if stop was requested
    if (status == 3)
        fprintf('\nFalsification stopped after max number of simulations.\n')       
        fprintf(['Final robustness: ' num2str(val_opt)]);
        return;    
    end
   
end

fprintf(['\nMax number of falsification restarts reached.\n']);
fprintf(['Final robustness: ' num2str(val_opt)]);
fprintf('\n');

end
