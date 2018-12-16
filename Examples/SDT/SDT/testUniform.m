close all;
clear;
clc;
warning('OFF', 'ALL')

addpath('./Functions')
addpath('./ExampleBand-Pass')
addpath('./ExampleLow-Pass')
addpath('./periodic.signals')
addpath('/Users/thaodang/breach-dev')



InitBreach
BP2param

model_name = 'BP2_in';

%rng(15000,'twister');
BrSD = BreachSimulinkSystem(model_name);
%, 'all', [], {}, [], 'Verbose',0,'SimInModelsDataFolder', true);


for sigId = 1:51
    if (sigId < 10)
        sigfilename = strcat('s_0', num2str(sigId));
    else
        sigfilename = strcat('s_', num2str(sigId));
    end
    
    In1 = load(sigfilename, '-ascii');
    scaling=1e-8;
    In1(:,1) = scaling*In1(:,1); 
    time = In1(:,1);
    time= time-time(1);
    sg_in = from_workspace_signal_gen({'In1'});
    
    BrSD_temp=BrSD.copy();
    BrSD_temp.SetInputGen({sg_in});
    BrSD_temp.Sim(time);
    
    
    %BrSD.SetTime([0 sim_time]);
    phi = STL_Formula('notsaturation', 'alw_[0, sim_time] (Out1[t] < 2)');
    phi = set_params(phi, {'sim_time'}, time(end));
    
    BrSD_temp.Sim(time);
    %input('Press ENTER to continue');
    %new_fig = gcf;
    BrSD_temp.PlotSignals({'In1', 'Out1'});
    
    
    rob = BrSD_temp.CheckSpec(phi)
    sigId
    if rob<0
        disp('Falsified in StatFalsify');
        
        %BrSD_temp.PlotSignals({'In1', 'OutSat'});
        falsified = true;
        sigId
        break;
    end
    
end


