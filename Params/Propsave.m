function Propsave(Sys, varargin)
%  
% PROPSAVE save properties set in the default properties set file of a system
%
% Synopsis:  Propsave(Sys, 'phi1', 'phi2', ..., )
%
%


  
  filename = [Sys.Dir filesep SysName(Sys) '_properties.mat'];

  for i = 1:numel(varargin)
    P___ = varargin{i};
    
    if (~isstr(P___))
      error(['Arguments should be a system structure and strings naming ' ...
             'parameter sets' ]);
    end
    
    eval([P___ '= evalin(''base'',''' P___ ''');']);
    if(exist(filename,'file')==0)
      save(filename, P___);
    else
      save(filename, '-append', P___);
    end
  end

