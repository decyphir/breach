%% Demonstrate signal generators reading files 

%% Create example files
%
BrDemo.InitAFC
Br0 = BrAFC.copy();
BrFromFiles=BrAFC.copy();

%% 
% Creates 3 simulations using the default input generator 
Br0.SetParamRanges('Pedal_Angle_base_value', [0 61.1]);
Br0.GridSample(3);
Br0.Sim()

%% Create files 
%  Files should be mat files with with a 'Time' array and a 'Pedal_Angle'
%  and 'Engine_Speed' arrays of the same dimensions 
%

Time = Br0.GetTime();
Pedal_Angle_values = Br0.GetSignalValues('Pedal_Angle'); % this has three individual traces 
Engine_Speed_values =  Br0.GetSignalValues('Engine_Speed'); % this has three individual traces 
for ifile = 1:3
    Pedal_Angle = Pedal_Angle_values{ifile}(1,:)'; 
    Engine_Speed = Engine_Speed_values{ifile}(1,:)'; 
    filename = [ 'AFC_input_data_' num2str(ifile) '.mat'];
    save(filename , 'Time', 'Pedal_Angle', 'Engine_Speed');
    filenames{ifile} = filename; 
end

%% Create signal generator reading files
% 
sg = from_file_signal_gen({'Pedal_Angle', 'Engine_Speed'}, filenames); 
BrFromFiles.SetInputGen(sg);

%%
% BrFromFile now reads inputs from files. It has one parameter named
% Pedal_Angle_file_idx which controls which file is read. 

BrFromFiles.SetParam('Pedal_Angle_file_idx', 2); % Computes using file 2
BrFromFiles.Sim(); 
BrFromFiles.PlotSignals({'Pedal_Angle', 'AF'});
 
 %%
 %
 BrFromFiles.ResetSimulations();
 BrFromFiles.SetParam('Pedal_Angle_file_idx',[1 3]); % computes using file 1 and 3
 BrFromFiles.Sim(); 
 figure; BrFromFiles.PlotSignals({'Pedal_Angle', 'AF'});