%% !!The model requires Matlab R2015a for some reason!!

%% Set thing up
clear;
InitBreach;

%% Solver config
my_solver = 'global_nelder_mead';
%my_solver = 'fmincon';
%my_solver = 'cmaes';
%my_solver = 'ga';

%% Configure and load the model
fuel_inj_tol = 1.0; 
MAF_sensor_tol = 1.0;
AF_sensor_tol = 1.0; 
pump_tol = 1.;
kappa_tol=1; 
tau_ww_tol=1;
fault_time=50;
kp = 0.04;
ki = 0.14;
mdl = 'AbstractFuelControl';

warning('off', 'Simulink:LoadSave:EncodingMismatch') % Not sure why, we just don't like warnings...

BrAFC = BreachSimulinkSystem(mdl, 'all', [], {}, [], 'Verbose',0,'SimInModelsDataFolder', false); 

%% Load specifications
STL_ReadFile('AFC_example_spec.stl');

%% Define a piecewise constant input signal generator
pedal_angle_gen = fixed_cp_signal_gen({'Pedal_Angle'}, ... % signal name
                                       3,...                % number of control points
                                      {'previous'});       % interpolation method 
                                  
engine_gen = constant_signal_gen({'Engine_Speed'}); 

InputGen = BreachSignalGen({pedal_angle_gen, engine_gen});   
                                  
InputGen.SetParam({'Pedal_Angle_u0','Pedal_Angle_u1','Pedal_Angle_u2','Engine_Speed_u0'},...
                   [30 30 30 1000]);
               
BrAFC.SetInputGen(InputGen);

BrAFC.SetParamRanges({'Pedal_Angle_u0', 'Pedal_Angle_u1','Pedal_Angle_u2','Engine_Speed_u0'}, [10 60; 10 60; 10 60; 900 1100 ]);

%% Falsify with classic robustness
display('---------------------------------- Classic Robustness -------------------------------------');
AFC_Falsify = BrAFC.copy();

AFC_Falsify = FalsificationProblem(AFC_Falsify, Overshoot_req);
AFC_Falsify.setup_solver(my_solver);
tic
AFC_Falsify.solve();
toc

%% Falsify without any guidance (use Boolean satisfaction value)
display('---------------------------------- Random -------------------------------------');
AFC_Falsify_Rand = BrAFC.copy();

AFC_Falsify_Rand = FalsificationProblem(AFC_Falsify_Rand, Overshoot_req);
AFC_Falsify_Rand.setup_solver(my_solver);
AFC_Falsify_Rand.set_IO_robustness_mode('random');
tic
AFC_Falsify_Rand.solve();
toc

%% Falsify with input robustness
display('---------------------------------- Input Robustness -------------------------------------');
AFC_Falsify_I = BrAFC.copy();

AFC_Falsify_I = FalsificationProblem(AFC_Falsify_I, Overshoot_req);
AFC_Falsify_I.setup_solver(my_solver);
AFC_Falsify_I.set_IO_robustness_mode('in',10^10); % fmincon does not like
                                                  % Inf, so we replace it
                                                  % by 10^10
tic
AFC_Falsify_I.solve();
toc

%% Falsify with output robustness
display('---------------------------------- Output Robustness -------------------------------------');
AFC_Falsify_O = BrAFC.copy();

AFC_Falsify_O = FalsificationProblem(AFC_Falsify_O, Overshoot_req);
AFC_Falsify_O.setup_solver(my_solver);
AFC_Falsify_O.set_IO_robustness_mode('out',10^10); % fmincon does not like
                                                   % Inf, so we replace it
                                                   % by 10^10
tic
AFC_Falsify_O.solve();
toc

%% State and Solve the same Falsification Problem, but with Combined IO Robustness
display('---------------------------------- Combined IO Robustness -------------------------------------');
AFC_Falsify_Rio = BrAFC.copy();

AFC_Falsify_Rio = FalsificationProblem(AFC_Falsify_Rio, Overshoot_req);
AFC_Falsify_Rio.setup_solver(my_solver);
AFC_Falsify_Rio.set_IO_robustness_mode('combined'); 
tic
AFC_Falsify_Rio.solve();
toc

%% State and Solve the same Falsification Problem, but with Constrained IO Robustness
display('---------------------------------- Constrained IO Robustness -------------------------------------');
AFC_Falsify_Cio = BrAFC.copy();

AFC_Falsify_Cio = FalsificationProblem(AFC_Falsify_Cio, Overshoot_req);
AFC_Falsify_Cio.setup_solver(my_solver);
AFC_Falsify_Cio.set_IO_robustness_mode('constrained',10^10); % fmincon does not like
                                                  % Inf, so we replace it
                                                  % by 10^10
tic
AFC_Falsify_Cio.solve();
toc

%% Plotting

% Options for 2D plot of parameter space
p1 = 2; % index of parameter 1
p2 = 3; % index of parameter 2

% Standard robustness and random falsification
figure;
subplot(2,2,1);
plot(AFC_Falsify.obj_log);
title('Standard cost function')
subplot(2,2,2);
plot(AFC_Falsify_Rand.obj_log);
title('Flat cost function')
subplot(2,2,3);
scatter(AFC_Falsify.X_log(p1,:),AFC_Falsify.X_log(p2,:),linspace(1,100,length(AFC_Falsify.X_log)));
title('Standard exploration of Pedal\_Angle');
subplot(2,2,4);
scatter(AFC_Falsify_Rand.X_log(p1,:),AFC_Falsify_Rand.X_log(p2,:),linspace(1,100,length(AFC_Falsify_Rand.X_log)))
title('Flat-earth exploration of Pedal\_Angle');

% Input and output robsutness
figure;
subplot(2,2,1);
plot(AFC_Falsify_I.obj_log);
title('Input cost function')
subplot(2,2,2);
plot(AFC_Falsify_O.obj_log);
title('Output cost function')
subplot(2,2,3);
scatter(AFC_Falsify_I.X_log(p1,:),AFC_Falsify_I.X_log(p2,:),linspace(1,100,length(AFC_Falsify_I.X_log)))
title('Input-driven exploration of Pedal\_Angle');
subplot(2,2,4);
scatter(AFC_Falsify_O.X_log(p1,:),AFC_Falsify_O.X_log(p2,:),linspace(1,100,length(AFC_Falsify_O.X_log)))
title('Output-driven exploration of Pedal\_Angle');

% Combined and constrained IO
figure;
subplot(2,2,1);
plot(AFC_Falsify_Rio.obj_log);
title('Combined IO cost function')
subplot(2,2,2);
plot(AFC_Falsify_Cio.obj_log);
title('Constrained IO cost function')
subplot(2,2,3);
scatter(AFC_Falsify_Rio.X_log(p1,:),AFC_Falsify_Rio.X_log(p2,:),linspace(1,100,length(AFC_Falsify_Rio.X_log)))
title('Combined IO-driven exploration of Pedal\_Angle');
subplot(2,2,4);
scatter(AFC_Falsify_Cio.X_log(p1,:),AFC_Falsify_Cio.X_log(p2,:),linspace(1,100,length(AFC_Falsify_Cio.X_log)))
title('Constrained IO-driven exploration of Pedal\_Angle');


