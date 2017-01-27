InitBreach;
%% load formula
formulas = STL_ReadFile('spec.stl');
phi_template = phi_100;
  
%% Create system and input strategy 
mdl = 'Autotrans_shift';
Br = BreachSimulinkSystem('Autotrans_shift');
Br.SetInputGen('VarStep2');
Br.SetTime(0:.01:50);

% Accelerate for 20 s at 100%
Br.SetParam('throttle_dt0', 20);
Br.SetParam('throttle_u0', 100);
Br.SetParam('brake_dt0', 20);
Br.SetParam('brake_u0', 0);
     
% Brake for 20s after
Br.SetParam('throttle_dt1', 20);
Br.SetParam('throttle_u1',   0);
Br.SetParam('brake_dt0', 20);
Br.SetParam('brake_u0', 325);

%% Property parameters 
prop_opt.params = {'t1', 'rpm_min'};
prop_opt.monotony = [-1 -1];

prop_opt.order = [2 1];
prop_opt.ranges = [0 50; ...
                       0 6000; ... 
                       ];
%% Property parameters 
prop_params.names = { 't1', 'rpm_min'};
prop_params.ranges = [0 50; ...
                       0 6000; ... 
                       ];
  
%% System and falsification parameters
input_params.names = {'throttle_dt0', 'throttle_u0','brake_dt0','brake_u1'};
input_params.ranges = [ 0 20   ;  ...  
                      0 100  ;  ...
                      0 20   ;
                      0 325 ;
                       ];
 
mine_phi_100 = ReqMiningProblem(Br, phi_template, input_params, prop_params);
mine_phi_100.synth_pb.solver_options.monotony= [-1 -1]; 
mine_phi_100.solve();