function Sys = CreateExternSystem(name, vars, params, p0, simfn)
% CREATEEXTERNSYSTEM creates a system with or without custom external simulator     
%
% Synopsis: Sys = CreateExternSystem(name, vars, params, p0, simfn)
%  
% Creates a system structure Sys to be used with Breach and a potential external
% simulator. Note: for Simulink, use CreateSimulinkSystem.  simfn must be a
% function with prototype  
%   function  [t X] = simfn(Sys, tspan, pts)  
% 
% Example:   
%
%   name = 'blackbox'     %  some name for the system
%
%   signals = {'s0','s1'};  % variables for signals values
%
%   params = {'p0','p1','p2'};  % system's parameters 
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
              vars = {vars{:} ['x' num2str(i)]};
          end
      end
      
      if isscalar(params)
          nparams = params;
          params = {};
          
          for i=1:nparams
              params = {params{:} ['x' num2str(i)]};
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
  else
    p0 = reshape(p0, numel(p0),1);
    if (numel(p0)== Sys.DimP)
      Sys.p = p0;
    else
      Sys.p = [Sys.x0'; p0];
    end
  end
  
  Sys.type = 'Extern';
  if exist('simfn','var')
      Sys.sim = simfn;
  end
  Sys.tspan= 0:.1:1;
  Sys.Dir = pwd;
