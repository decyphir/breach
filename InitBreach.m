% This script initializes Breach, in particular adding paths to Breach directories

cdr = pwd;
dr = which('InitBreach');
dr =  dr(1:end-13);

cd(dr);

if (exist('BreachGlobOpt.mat'))  
  load BreachGlobOpt;
else
  BreachGlobOpt.breach_dir = dr;  
  save('BreachGlobOpt.mat', 'BreachGlobOpt')
end

addpath(dr);
addpath( [dr filesep 'Core']);
addpath( [dr filesep 'Core' filesep 'm_src']);
addpath( [dr filesep 'Params']);
addpath( [dr filesep 'Params' filesep 'm_src']);
addpath( [dr filesep 'Plots']);
addpath( [dr filesep 'Plots' filesep 'm_src']);
addpath( [dr filesep 'Toolboxes' filesep 'optimize']);
addpath( [dr filesep 'Toolboxes' filesep 'sundials' filesep 'sundialsTB' ]);
addpath( [dr filesep 'Toolboxes' filesep 'sundials' filesep 'sundialsTB' filesep 'cvodes']);

%addpath(genpath(BreachGlobOpt.breach_dir));

%% Init BreachGlobOpt options and fourre-tout global variable
BreachGlobOptTmp = BreachGlobOpt;
clear BreachGlobOpt;
global BreachGlobOpt;
BreachGlobOpt = BreachGlobOptTmp; 
clear BreachGlobOptTmp;

BreachGlobOpt.RobustSemantics = 0 ; % 0 by default, -1 is for left time robustness, +1 for right, inf for sum ?

cd(cdr);
