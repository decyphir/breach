function Psave(Sys, varargin)
%  
% PSAVE save a parameter set in the default param set file of a system
%
% Synopsis:  Psave(Sys, 'P1', 'P2', ..., )
%
% Ex: 
%
%  P = CreateSampling(Sys);
%  Ph = QuasiRefine(P,100);
%  P2 = Refine(P,2);
%  Psave(Sys, 'P', 'Ph', 'P2' );
%


  
  filename = [Sys.Dir filesep SysName(Sys) '_param_sets.mat'];

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

