% This script initializes Breach, in particular adding paths to Breach directories

cdr = pwd;
dr = which('InitBreach');
dr =  dr(1:end-13);

cd(dr);

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

if (exist('BreachGlobOpt.mat'))  
  load BreachGlobOpt;
  
  % Convert BreachGlobOpt into global
  BreachGlobOptTmp = BreachGlobOpt;
  clear BreachGlobOpt;
  global BreachGlobOpt;
  BreachGlobOpt = BreachGlobOptTmp; 
  clear BreachGlobOptTmp;
  BreachGlobOpt.RobustSemantics = 0 ; % 0 by default, -1 is for left time robustness, +1 for right, inf for sum ?

else 
  if (~exist('BreachGlobOpt','var'))
    global BreachGlobOpt;
    BreachGlobOpt.breach_dir = dr;  
    BreachGlobOpt.RobustSemantics = 0;
  end
  
  if (~isfield(BreachGlobOpt,'RobustSemantics'))        
    BreachGlobOpt.RobustSemantics = 0;  
  end
end

cd(cdr);
