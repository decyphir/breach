function Psave(Sys, varargin)
%PSAVE saves one or many parameter set in the default param set file of a
% system. There is two ways to use Psave. The first one consist to provide
% in the arguments the names of the parameter sets to save. In this case,
% the parameter sets must be defined in the workspace (ie: if they are
% defined inside a function, they are not in the workspace). The other way
% to use this function consists to provide in the second argument the name
% of the parameter set, and in the third, the parameter set itself. In this
% case, the parameter set do not have to be defined in the workspace.
% 
% Synopsis:  Psave(Sys, 'P1', 'P2', ...)
%       OR:  Psave(Sys, 'Pname', Pvar)
% 
% Inputs:
%  - Sys     : the system
%  - 'P1'    : strings describing the name of the parameter set to save
%  - 'Pname' : string describing the name with which the parameter set is
%              saved
%  - Pvar    : parameter set to save
% 
% Example (Lorentz84):
%   P = CreateParamSet(Sys);
%   Ph = QuasiRefine(P,100);
%   P2 = Refine(P,2);
%   Psave(Sys, 'P', 'Ph', 'P2');
% 
% Other example:
%   function myFunc(Sys)
%       P = CreateParamSet(Sys);
%       Ph = QuasiRefine(P,100);
%       Psave(Sys, 'Ph', Ph);
%   end
% 
%See also Propsave
%

filename = [Sys.Dir filesep SysName(Sys) '_param_sets.mat'];

% Case 1: we provide the name and the parameters set to save
if(numel(varargin)==2)
    if(ischar(varargin{1}) && isstruct(varargin{2}))
        eval([varargin{1},'=varargin{2};']);
        if(exist(filename,'file')==0)
            save(filename, varargin{1});
        else
            save(filename, '-append', varargin{1});
        end
        return ;
    end
end

% Case 2: the parameters sets to save are defined in the environment, and
%           only their names are provided
for ii = 1:numel(varargin)
    P___ = varargin{ii};
    
    if ~ischar(P___)
        error('Psave:parameters',['Arguments should be a system ' ...
            'structure and strings naming parameter sets' ]);
    end
    
    eval([P___ '= evalin(''base'',''' P___ ''');']);
    if(exist(filename,'file')==0)
        save(filename, P___);
    else
        save(filename, '-append', P___);
    end
end

end
