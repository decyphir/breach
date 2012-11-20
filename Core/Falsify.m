function Pr = Falsify (Sys, phi, falsif_opt, params_prop)
%
%  FALSIFY Tries to falsify property phi 
%
%  Synopsis:    Pr = Falsify(Sys, phi,  falsif_opt [, param_prop])
%
%  Example:
%  
%  mdl = 'Autotrans_shift';
%  Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'UniStep1');
%  Sys.tspan = 0:.01:50;
%  
%  falsif_prop.params_u = {'throttle_u0'};
%  falsif_opt.range = [0 100];
%
%  falsif_opt.nb_init = 10;
%  falsif_opt.nb_iter = 100;
%  
%  Pr = Falsify(Sys, phi1, falsif_opt);
%  
  
%% Create system and input strategy
  
  params =falsif_opt.params;
  ranges = falsif_opt.ranges;
  Pu = CreateParamSet(Sys,params, ranges);
 
  if exist('params_prop','var')
    Pu = SetParam(Pu, params_prop.names, ...
                      params_prop.values);
  end
      
  nb_init = falsif_opt.nb_init;
  nb_iter = falsif_opt.nb_iter;
      
  %% Arguments for falsification
  
  opt.OptimType = 'Min';
  opt.MaxIter = nb_iter;
  opt.StopWhenFoundInit = 1;
  opt.StopWhenFound = 1;
  opt.ubound = ranges(:,2);
  opt.lbound = ranges(:,1);
 
  % Get a set of initial search points     
  Pr_new = QuasiRefine(Pu, nb_init);  
    
  % Falsification step 
  [val_opt Propt] = SOptimProp(Sys, Pr_new, phi, opt);
  ifalse = find(val_opt<0);
    
  if ~isempty(ifalse)
    Pfals = Sselect(Propt, ifalse);        
    fprintf('  ---- Falsified !');    
    Pf = SConcat(Pr_new,Pfals);
  else 
    fprintf(['\nStopped. Final robustness: ' num2str(val_opt)]);  
    fprintf('\n');
    return
  end
  fprintf('\n');
  
 