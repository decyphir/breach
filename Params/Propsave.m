function Propsave(Sys, varargin)
%PROPSAVE saves properties set in the default properties set file of a system
%
% Synopsis:  Propsave(Sys, 'phi1', 'phi2', ..., )
%
%See also Psave
%

global BreachGlobOpt
filename = [Sys.Dir filesep SysName(Sys) '_properties.mat'];

if(nargin==2)
    if isstruct(varargin{1})
        struct2save= varargin{1};
    end
    
    if ischar(varargin{1})
        struct2save = struct(varargin{1}, BreachGlobOpt.STLDB(varargin{1}));
    end
elseif nargin>=2
    for i = 1:numel(varargin)
        if ~ischar(varargin{i})
            error('Propsave:argumentError',['Arguments should be a system structure and strings naming ' ...
                'properties or the name of a structure with properties.']);
        end
        phi_id = varargin{i};
        struct2save.(phi_id) = BreachGlobOpt.STLDB(varargin{i});
    end
else
    error('Propsave:argumentError',['Arguments should be a system structure and strings naming ' ...
        'properties or the name of a structure with properties.']);
end

if ~exist(filename,'file')
    save(filename, '-struct', 'struct2save');
else
    save(filename, '-struct', 'struct2save', '-append');
end

end
