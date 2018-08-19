% Init Model with two inputs 

mdl = 'Autotrans_shift';
BrAutotrans_nominal = BreachSimulinkSystem(mdl, 'all', [], {}, [], 'Verbose',0,'SimInModelsDataFolder', true); 

BrAutotrans_nominal.SetTime(0:.01:40); % default simulation time
BrAutotrans_nominal.SetInputGen('VarStep2') 

%% Set input values (other than 0) 
% Accelerate for 20 s at 100%
BrAutotrans_nominal.SetParam( 'throttle_dt0', 20);
BrAutotrans_nominal.SetParam( 'throttle_u0', 100);
BrAutotrans_nominal.SetParam( 'brake_u0', 0);
BrAutotrans_nominal.SetParam( 'brake_dt0', 20);
     
% Brake ever after
BrAutotrans_nominal.SetParam( 'throttle_u1', 0);
BrAutotrans_nominal.SetParam( 'brake_u1', 325);

%% 3 Ranges
BrAutotrans_3ranges = BrAutotrans_nominal.copy();

% Acceleration
BrAutotrans_3ranges.SetParamRanges( 'throttle_u0', [0 100]);
BrAutotrans_3ranges.SetParam( 'throttle_dt0', 20);
BrAutotrans_3ranges.SetParamRanges( 'throttle_u1', [0 100]);
     
% Braking
BrAutotrans_3ranges.SetParam('brake_u0', 0);
BrAutotrans_3ranges.SetParam('brake_dt0', 20);
BrAutotrans_3ranges.SetParamRanges( 'brake_u1', [0 325]);

%% 6 ranges
BrAutotrans_6ranges = BrAutotrans_nominal.copy();

% Acceleration
BrAutotrans_6ranges.SetParamRanges( 'throttle_u0', [0 100]);
BrAutotrans_6ranges.SetParamRanges( 'throttle_dt0', [0 40]);
BrAutotrans_6ranges.SetParamRanges( 'throttle_u1', [0 100]);
     
% Braking
BrAutotrans_6ranges.SetParamRanges( 'brake_u0', [0 325]);
BrAutotrans_6ranges.SetParamRanges( 'brake_dt0', [0  40]);
BrAutotrans_6ranges.SetParamRanges( 'brake_u1', [0 325]);

%% 
BrAutotrans_varstep5 =BreachSimulinkSystem(mdl, 'all', [], {}, [], 'Verbose',0,'SimInModelsDataFolder', true); 
sg = var_step_signal_gen({'throttle', 'brake'}, 5);
BrAutotrans_varstep5.SetInputGen(sg);
            
% We assign ranges for duration and amplitude of each input:
BrAutotrans_varstep5.SetParamRanges({'dt_u0', 'dt_u1', 'dt_u2', 'dt_u3'}, ...
                  [.1 10  ;  .1 10;    0.1 10;    0.1 10]);
BrAutotrans_varstep5.SetParamRanges({'throttle_u0','brake_u1', 'throttle_u2', 'brake_u3'}, ... 
                  [0 100;        0 325;      0 100;         0 325]);
