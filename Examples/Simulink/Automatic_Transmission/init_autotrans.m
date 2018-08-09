% Init Model with two inputs 

mdl = 'Autotrans_shift';
Br = BreachSimulinkSystem(mdl);
Br.SetTime(0:.01:40); % default simulation time
Br.SetInputGen('VarStep2') 

Br.PrintAll()

% Set input values (other than 0) 
% Accelerate for 20 s at 100%
Br.SetParam( 'throttle_dt0', 20);
Br.SetParam( 'throttle_u0', 100);
Br.SetParam( 'brake_dt0', 20);
Br.SetParam( 'brake_u0', 0);
     
% Brake for 20s after
Br.SetParam( 'throttle_u1', 0);
Br.SetParam( 'brake_dt0', 20);
Br.SetParam( 'brake_u0', 0);
Br.SetParam( 'brake_u1', 325);
