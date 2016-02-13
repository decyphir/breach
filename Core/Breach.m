function Breach(Sys, P, spec)
% BREACH  Launches Breach main GUI.
%
%   Usage: Breach([Sys, P, spec_file]) or Breach([spec_file, P])
%
%   Launches the main GUI of Breach (legacy function). If no system Sys is given as an
%   argument, it first checks whether a system is present in the current
%   directory (through the presence of CreateSystem.m file, and if it is not
%   the case, it proposes to create one. P is a parameter set name in the
%   workspace (string), spec_file is the name of a file with STL formulas.
%

% load Breach options
try
    BreachGlobOpt =  evalin('base', 'BreachGlobOpt');
catch
    evalin('base', 'InitBreach');
end

% Checks if a system is present in the current directory

dr = pwd;
indst = strfind(dr, filesep);
SysName = dr(indst(end)+1:end);

if (exist('Sys'))
    
    % if Sys is a string, try see if it's a Simulink model
    if ischar(Sys)
        try % is this a simulink model?
            Sys = CreateSimulinkSystem(Sys);
            Breach(Sys);
        catch
            try % is this an stl file
                [prop_names, props, signal_names] = STL_ReadFile(Sys);
                [~, Sys_name, ~] = fileparts(Sys);
                Sys = CreateExternSystem(Sys_name, signal_names, {},[]);
                Propsave(Sys, prop_names{:});
            catch
                error('Breach: system not recognized');
            end
        end
    end
    
    if ~isfield(Sys, 'type')
        Sys.type= 'Breach';
    end
    if strcmp(Sys.type, 'Simulink')
        SysName = Sys.mdl;
        load_system(Sys.mdl);
    end
    Sys.Dir= pwd;
    
    %% Extra arguments
    if nargin>=2 % a parameter set is given as input, load it
        if ~(isempty(P))
            Psave(Sys, P);
        end
    end
    
    if nargin>=3
        if ischar(spec) % 
            Propsave(Sys, spec);
        else % should be a cell
            Propsave(Sys, spec{:});        
        end
    end
    
    h = BreachGui('varargin',Sys,SysName);
    return
    
end

% either SysName.mat does not exist or not valid, now try CreateSystem
if (exist([pwd '/CreateSystem.m'])~=2)
    fprintf('No system found in the current directory.\nYou can set up one using the function NewSystem\nor CreateExternSystem if using another simulator\nor CreateSimulinkSystem for Simulink models \n');
    return;
else
    CreateSystem;
end

Sys.Dir = pwd;
save([SysName '.mat'],'Sys');

%% Extra arguments
if nargin>=2 % a parameter set is given as input, load it
    if ~(isempty(P))
        try
            Psave(Sys, P);
        catch
            warning(['Problem loading parameter set' P]);
        end
    end
end

if nargin>=3
    try
        prop_names = STL_ReadFile(spec);
        Propsave(Sys, prop_names{:});
    catch
        warning(['Problem loading properties from file ' spec]);
    end
end

BreachGui('varargin',Sys,SysName);
