clear;
InitBreach;

% AFC parameters
fuel_inj_tol = 1.0; 
MAF_sensor_tol = 1.0;
AF_sensor_tol = 1.0; 
pump_tol = 1.;
kappa_tol=1; 
tau_ww_tol=1;
fault_time=50;
kp = 0.04;
ki = 0.14;

warning('off', 'Simulink:LoadSave:EncodingMismatch')
mdl = 'AFC_Online';

max_rob = 20;
sig_names = 'AF,AFref';
phi_st = 'AF[t]>0';
phi_stop = 'AF[t]>0';

BrAFC_Online = BreachSimulinkSystem(mdl);

pedal_angle_gen = pulse_signal_gen({'Pedal_Angle'}); % Generate a pulse signal for pedal angle
engine_gen      = fixed_cp_signal_gen({'Engine_Speed'}, ... % signal name
                                       3,...                % number of control points
                                      {'spline'});       % interpolation method 
        
InputGen = BreachSignalGen({pedal_angle_gen, engine_gen});

InputGen.SetParam({'Engine_Speed_u0','Engine_Speed_u1','Engine_Speed_u2'},...
                        [1000 1100 500]);

InputGen.SetParam({'Pedal_Angle_base_value', 'Pedal_Angle_pulse_period', ...
                         'Pedal_Angle_pulse_amp','Pedal_Angle_pulse_width'}, ... 
                         [0 15 30 .5]);

BrAFC_Online.SetInputGen(InputGen);

% online block init
pedal_angle_gen = pulse_signal_gen({'Pedal_Angle'}); % Generate a pulse signal for pedal angle
engine_gen      = fixed_cp_signal_gen({'Engine_Speed'}, ... % signal name
                                       3,...                % number of control points
                                      {'spline'});       % interpolation method 
        
InputGen = BreachSignalGen({pedal_angle_gen, engine_gen});

InputGen.SetParam({'Engine_Speed_u0','Engine_Speed_u1','Engine_Speed_u2'},...
                        [1000 1100 500]);

InputGen.SetParam({'Pedal_Angle_base_value', 'Pedal_Angle_pulse_period', ...
                         'Pedal_Angle_pulse_amp','Pedal_Angle_pulse_width'}, ... 
                         [0 15 30 .5]);
   
BrAFC_Online.SetInputGen(InputGen);
