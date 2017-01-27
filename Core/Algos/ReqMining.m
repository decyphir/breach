function [p, rob, Pr] = ReqMining(Sys, phi, falsif_opt, prop_opt, iter_max)
%REQMINING requirement mining algorithm
%
% Synopsis: [p, rob, Pfals] = Req_Mining(Sys, template_phi, falsif_opt, prop_opt, iter_max)
%
%  Example:
%  
%  mdl = 'Autotrans_shift';
%  Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'UniStep1');
%  Sys.tspan = 0:.01:50;
%  
%  prop_opt.params = {'vmax', 'rpm_max'};
%  prop_opt.monotony   = [1 1];         % either +1 or -1
%                                       % the order should be consistent
%                                       % with field params
%            
%  prop_opt.p_interval = [0 200   ;...  % for vmax
%                         0 6000 ];     % for rmp_max
%  prop_opt.p_tol      = [1 1];
%  
%  falsif_opt.params = {'throttle_u0'};
%  falsif_opt.ranges = [0 100];
%
%  falsif_opt.iter_max = 100;
%  falsif_opt.nb_init = 10;
%  falsif_opt.nb_iter = 100;
%  
%  [p, rob,Pr] = ReqMining(Sys, phi1, falsif_opt, prop_opt);
%  
%  See also GetPropParamBin 
    
  Pprop = CreateParamSet(Sys, prop_opt.params);
 
  if (~exist('iter_max','var'))
      iter_max = falsif_opt.nb_iter;
  end 
        
  %% First trajectory
  
  i = 0;
  Pr = ComputeTraj(Sys, Pprop, Sys.tspan); % Pr stores the falsifying trajectories
  
  % Mine parameters for first simulation
  [p, rob] = GetPropParamBin(Sys, phi, Pr, prop_opt);
  
  prop_opt.values = p;
  i = 1;
  Pr = SetParam(Pr, prop_opt.params, prop_opt.values');
  
  iter = 1;
  fprintf('\n');
  
  %% Main loop
  while iter<=iter_max
    
    % update display
    clc;    
    fprintf('Iter %d/%d\n\n', iter, iter_max);
    fprintf('Attempts falsifying PSTL formula:\n\n');
    fprintf('%s\n\n', disp(phi)); 
    fprintf('with:\n\n');
    for i = 1:numel(prop_opt.params)
      fprintf('%s= %g\n', prop_opt.params{i}, p(i))
    end
    fprintf('\n');

    % Falsification step 
    [Propt, rob] = Falsify(Sys, phi, falsif_opt, prop_opt);
    
    ifalse = find(rob<0);
    
    if ~isempty(ifalse)
        Pfals = Sselect(Propt, ifalse);
        fprintf('\nGetting new parameters for formula:\n');
        Pr = SConcat(Pr,Pfals);
        % Parameter inferrence step
        [p, rob] = GetPropParamBin(Sys, phi, Pr,  prop_opt);
        prop_opt.values = p;
    else
        fprintf(['\np:' num2str(p) '\n']) ;  
        return;
    end
    iter=iter+1;
  end
  fprintf('\nReq mining stopped after max number of iterations.\n');
  fprintf(['\np: ' num2str(p) '\n' ]);
  
 
