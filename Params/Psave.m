function Psave(Sys, varargin)
%  
% PSAVE save a parameter set in the default param set file of a system
%
% Synopsis:  Psave(Sys, P1, P2, ..., )
%
% Ex: 
%
%  P = CreateSampling(Sys);
%  Ph = HaltonRefine(P,100);
%  P2 = Refine(P,2);
%  Psave(Sys, 'P', 'Ph', 'P2' );
%

  filename = [Sys.Dir filesep SysName(Sys) '_param_sets.mat'];

  for i = 1:numel(varargin)
    P___ = varargin{i};
    
    eval([P___ '= evalin(''base'',''' P___ ''');']);
    save(filename, '-append', P___);
  end

