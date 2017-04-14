%% Demonstrate signal generators reading files 

%% Create example files
%
BrDemo.InitAFC

%%
%
BrAFC.SetParamRanges('Pedal_Angle_base_value', [0 61.1]);

%% 
% 
BrAFC.GridSample(3);
BrAFC.Sim()

%% Create files
%
Pedal_Angle_values = BrAFC.GetSignalValues('Pedal_Angle');
Time = BrAFC.GetTime();

%%
%  Files should be mat files with variables being arrays of size N x 2, the
%  first column being time values and the second signal values.
for ifile = 1:3
    Pedal_Angle = [Time' Pedal_Angle_values{ifile}(1,:)'];  
    save( [ 'AFC_input_data_' num2str(ifile) '.mat'],  'Pedal_Angle');
end

%% Create signal generators
%
sg_pedal = from_file_signal_gen({'Pedal_Angle'}, 'AFC_input_data_*.mat',  'Pedal_Angle'); 
sg_engine = constant_signal_gen('Engine_Speed');
BrAFC.SetInputGen( {sg_pedal, sg_engine});

%%
% BrAFC now reads  pedal_angle input from files. 