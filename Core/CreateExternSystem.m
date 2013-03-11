function Sys = CreateExternSystem(name, vars, params, p0, simfn)
% CREATEEXTERNSYSTEM creates a system with or without custom external simulator     
%
% Synopsis: Sys = CreateExternSystem(name, vars, params, p0)
%  
% Creates a system structure Sys to be used with Breach and an external
% simulator. Note: for Simulink, use CreateSimulinkSystem.  
% 
% Example:   
%
%   signals = {'s0','s1'};  % variables for signals values
%
%   params = {'p0','p1','p2'};  % parameters related to signal or to be
%                               % used in temporal logic formulas
%
%   p0 = [0 0 0 ];    % default for parameters p0,p1 and p2
%
%   Sys = CreateExternSystem(name, signals ,params, p0); % creates the Sys structure
%
  
  Sys.name = name;

  if nargin==3
      if isscalar(vars)
          nvars = vars;
          vars = {};
          for i=1:nvars
              vars = {vars{:} ['y' num2str(i)]};
          end
      end
      
      if isscalar(params)
          nparams = params;
          params = {};
          
          for i=1:nparams
              params = {params{:} ['y' num2str(i)]};
          
          end          
      end
  end
  
  
  Sys.DimX = numel(vars);
  Sys.DimU =0; 
  Sys.DimP = numel(vars)+numel(params); 
  Sys.ParamList = {vars{:} params{:}};
  Sys.x0 = zeros(1,numel(vars));  
  
  if (~exist('p0'))
      Sys.p = zeros(1, Sys.DimP);
  elseif (numel(p0)==Sys.DimP)
      Sys.p = p0;
  else
      Sys.p = [Sys.x0 p0];
  end
  
  if exist('simfn','var')
      Sys.type = 'Extern';
      Sys.sim = simfn;
  else
      Sys.type = 'traces';
  end
  
  Sys.Dir = pwd;