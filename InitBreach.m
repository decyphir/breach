% This script initializes Breach, in particular adding paths to Breach directories
disp('Initializing Breach...'); 
disp('%-------------------------------------------------------------------%')
disp('% Important!                                                        %')
disp('% The prefix ''QMITL'' has (finally) been replaced by ''STL''.      %');
disp('% Please update your code accordingly if needed.                    %'); 
disp('%  - Find/replace QMITL with STL in all m-files will do the trick.) %');
disp('%  - Use RecoverFormula to convert previously saved QMITL_Formula   %')
disp('%    Object into STL_Formula.                                       %')
disp('%-------------------------------------------------------------------%')

id = 'MATLAB:dispatcher:nameConflict';
warning('off',id);

cdr_ = pwd;
dr_ = which('InitBreach');
dr_ =  dr_(1:end-13);

cd(dr_);

addpath(dr_);
addpath( [dr_ filesep 'Core']);
addpath( [dr_ filesep 'Core' filesep 'm_src']);
addpath( [dr_ filesep 'Params']);
addpath( [dr_ filesep 'Params' filesep 'm_src']);
addpath( [dr_ filesep 'Params' filesep 'm_src' filesep 'sobolqr']);
addpath( [dr_ filesep 'Params' filesep 'm_src' filesep 'niederreiter2']);
addpath( [dr_ filesep 'Plots']);
addpath( [dr_ filesep 'Plots' filesep 'm_src']);
addpath( [dr_ filesep 'Toolboxes' filesep 'optimize']);
addpath( [dr_ filesep 'Toolboxes' filesep 'sundials' filesep 'sundialsTB' ]);
addpath( [dr_ filesep 'Toolboxes' filesep 'sundials' filesep 'sundialsTB' filesep 'cvodes']);

%% Init BreachGlobOpt options and fourre-tout global variable

if exist('BreachGlobOpt.mat')
    load BreachGlobOpt;
    
    % Convert BreachGlobOpt into global
    BreachGlobOptTmp = BreachGlobOpt;
    clear BreachGlobOpt;
    global BreachGlobOpt;
    BreachGlobOpt = BreachGlobOptTmp;
    clear BreachGlobOptTmp;
    BreachGlobOpt.RobustSemantics = 0 ; % 0 by default, -1 is for left time robustness, +1 for right, inf for sum ?
    
else
    
    if ~exist('BreachGlobOpt','var')
        global BreachGlobOpt;
        BreachGlobOpt.breach_dir = dr_;
    end
    
    if ~isfield(BreachGlobOpt,'RobustSemantics')
        BreachGlobOpt.RobustSemantics = 0;
    end
    
end
cd(cdr_);
clear cdr_ dr_;

warning('on',id);
