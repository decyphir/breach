clear;

BrDemo.InitAFCparams;

mdl = 'AbstractFuelControl';

BrAFC = BreachSimulinkSystem(mdl, 'all', [], {}, [], 'Verbose',0,'SimInModelsDataFolder', true); 

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

BrAFC.SetInputGen(InputGen);
 