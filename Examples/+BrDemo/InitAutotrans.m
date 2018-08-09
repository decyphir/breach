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
