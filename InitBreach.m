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
addpath( [dr filesep 'Params']);
addpath( [dr filesep 'Plots']);
addpath( [dr filesep 'Toolboxes' filesep 'optimize']);
addpath( [dr filesep 'Toolboxes' filesep 'sundials' filesep 'sundialsTB' ]);
addpath( [dr filesep 'Toolboxes' filesep 'sundials' filesep 'sundialsTB' filesep 'cvodes']);

%addpath(genpath(BreachGlobOpt.breach_dir));

cd(cdr);