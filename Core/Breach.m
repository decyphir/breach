function gui = Breach(Sys, P, spec)
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

if nargin==1
    if (isa(Sys, 'BreachSystem'))
        BrSys = Sys;
        Sys = BrSys.Sys;
        P = BrSys.P;
    end
end

% Checks if a system is present in the current directory
dr = pwd;
indst = strfind(dr, filesep);
SysName = dr(indst(end)+1:end);

% if no BreachSystem is given, create one.


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
    
    % Extra arguments
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
    
else
   
    % either SysName.mat does not exist or not valid, now try CreateSystem
    if (exist([pwd '/CreateSystem.m'])~=2)
        error('No system found as argument or in the current directory.');
    else
        CreateSystem; % old stuff 
    end
    
    Sys.Dir = pwd;
    save([SysName '.mat'],'Sys');
    
    % Extra arguments
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
    
end


if ~exist('BrSys', 'var')
    BrSys = BreachSystem(Sys);
end
    

gui = BreachGui('varargin',Sys,SysName, BrSys);
