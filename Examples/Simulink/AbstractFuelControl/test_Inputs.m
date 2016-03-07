%% pulse and constant interpolation
pedal_angle_gen = pulse_signal_gen({'Pedal_Angle'});
engine_gen      = fixed_cp_signal_gen({'Engine_Speed'}, ... % signal name
                                       3,...                % number of control points
                                      {'previous'}...       % interpolation method 
                                     );

InputGen = BreachSignalGen({pedal_angle_gen, engine_gen});
InputGen.PrintParams();

InputGen.SetParam({'Engine_Speed_u0','Engine_Speed_u1','Engine_Speed_u2'},...
                        [100 1100 900]);
InputGen.SetParam({'Pedal_Angle_base_value', 'Pedal_Angle_pulse_period', ...
                         'Pedal_Angle_pulse_amp','Pedal_Angle_pulse_width'}, ... 
                         [0 15 30 .5]);
InputGen.PrintParams();

InputGen.Sim([0 10]);
InputGen.PlotSignals();
title('pulse and constant interpolation')

%% linear interpolation
pedal_angle_gen = pulse_signal_gen({'Pedal_Angle'});
engine_gen      = fixed_cp_signal_gen({'Engine_Speed'}, ... % signal name
                                       3,...                % number of control points
                                      {'linear'}...         % interpolation method 
                                     );

InputGen2 = BreachSignalGen({pedal_angle_gen, engine_gen});
InputGen2.PrintParams();

InputGen2.SetParam({'Engine_Speed_u0','Engine_Speed_u1','Engine_Speed_u2'},...
                        [100 1100 900]);
InputGen2.SetParam({'Pedal_Angle_base_value', 'Pedal_Angle_pulse_period', ...
                         'Pedal_Angle_pulse_amp','Pedal_Angle_pulse_width'}, ... 
                         [0 15 30 .5]);
InputGen2.PrintParams();
InputGen2.Sim([0 10]);
InputGen2.PlotSignals();
title('step and linear interpolation')



%% test with system
BrAFC.SetInputGen(InputGen)

BrAFC.Sim(40);
BrAFC.PlotSignals({'Pedal_Angle', 'Engine_Speed','cyl_fuel', 'AF'});

%% Test my_signal_generator
% 
my_input_gen = my_signal_generator(2);        % creates an instance with lambda=2
MyInputGen = BreachSignalGen({my_input_gen}); % Makes it a Breach system 

BrAFC_MyGen = BrAFC.copy();
BrAFC_MyGen.SetInputGen(MyInputGen);          % Plug it to the Simulink model
BrAFC_MyGen.PrintParams();                    % Makes sure the new parameters are visible


BrAFC_MyGen.Sim(40);
BrAFC_MyGen.PlotSignals({'Pedal_Angle', 'Engine_Speed','cyl_fuel', 'AF'});
