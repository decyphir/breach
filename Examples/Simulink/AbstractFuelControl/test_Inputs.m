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
                         [0 5 30 .5]);
InputGen.PrintParams();

InputGen.Sim([0 10]);
figure;
set(gcf, 'Name','Pulse and constant interpolation');
InputGen.PlotSignals();

%% linear interpolation
engine_gen      = fixed_cp_signal_gen({'Engine_Speed'}, ... % signal name
                                       3,...                % number of control points
                                      {'linear'}...         % interpolation method 
                                     );
                                 
InputGen2 = BreachSignalGen({pedal_angle_gen, engine_gen});
InputGen2.PrintParams();

InputGen2.SetParam({'Pedal_Angle_base_value', 'Pedal_Angle_pulse_period', ...
                         'Pedal_Angle_pulse_amp','Pedal_Angle_pulse_width'}, ... 
                         [0 5 30 .5]);
InputGen2.SetParam({'Engine_Speed_u0','Engine_Speed_u1','Engine_Speed_u2'},...
                        [100 1100 900]);
InputGen2.PrintParams();
InputGen2.Sim([0 10]);
figure;
set(gcf, 'Name','Pulse and Linear Interpolation');
InputGen2.PlotSignals();

%% linear and spline interpolation with var steps 

input_gen = var_cp_signal_gen({'Engine_Speed', 'Pedal_Angle'}, ... % signal names
                                       [3 6],...                   % number of control points
                                       {'spline', 'linear'}...     % interpolation method
                                     );

InputGen3 = BreachSignalGen({input_gen});
InputGen3.PrintParams();

InputGen3.SetParam({'Pedal_Angle_u0','Pedal_Angle_dt0',...
                    'Pedal_Angle_u1','Pedal_Angle_dt1',...
                    'Pedal_Angle_u2','Pedal_Angle_dt2',...
                    'Pedal_Angle_u3','Pedal_Angle_dt3',...
                    'Pedal_Angle_u4','Pedal_Angle_dt4',...
                    'Pedal_Angle_u5'},...
                    [20 1 10 1 80 2 40 3 10 2.5 15]);

InputGen3.SetParam({'Engine_Speed_u0','Engine_Speed_dt0',...
                    'Engine_Speed_u1','Engine_Speed_dt1',...
                    'Engine_Speed_u2'}, ...
                    [1100 5 900 5 1000]);

InputGen3.PrintParams();
InputGen3.Sim([0 10]);
figure;
set(gcf, 'Name','Interpolation with variable steps');
InputGen3.PlotSignals();

%% test with system
BrDemo.InitAFC;

BrAFC.SetInputGen(InputGen)
BrAFC.Sim(40);
figure;
BrAFC.PlotSignals({'Pedal_Angle', 'Engine_Speed','cyl_fuel', 'AF'});

%% Test my_signal_generator
% 
my_input_gen = my_signal_generator(2);        % creates an instance with lambda=2
MyInputGen = BreachSignalGen({my_input_gen}); % Makes it a Breach system 

BrAFC_MyGen = BrAFC.copy();
BrAFC_MyGen.SetInputGen(MyInputGen);          % Plug it to the Simulink model
BrAFC_MyGen.PrintParams();                    % Makes sure the new parameters are visible

BrAFC_MyGen.Sim(40);
figure;
BrAFC_MyGen.PlotSignals({'Pedal_Angle', 'Engine_Speed','cyl_fuel', 'AF'});
