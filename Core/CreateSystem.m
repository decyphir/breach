function Sys = CreateSystem(vars, params, p0)
% CREATESYSTEM creates a dynamicless system      
%
% Synopsis: Sys = CreateSystem(vars, params, p0)
%  
% Example:   
%
%   vars = {'s0','s1'};  % variables for signals values
%
%   params = {'p0','p1','p2'};  % parameters related to signal or to be
%                               % used in temporal logic formulas
%
%   p0 = [ 0 0 ...      % default for initial values of signals 
%          0 0 0 ]';    % default for parameters p0,p1 and p2
%
%   Sys = CreateSystem(vars,params, p0); % creates the Sys structure
%
  
  
  Sys.DimX = numel(vars);
  Sys.DimU =0; 
  Sys.DimP = numel(vars)+numel(params); 
  Sys.ParamList = {vars{:} params{:}};
  Sys.x0 = zeros(1,numel(vars));  
  
  if (~exist('p0'))
      Sys.p = zeros(1, Sys.DimP);
  else
      Sys.p = p0;
  end
  
  Sys.x0(1:Sys.DimX)= p0(1:Sys.DimX);
  
  % dynamicless system
  Sys.type = 'traces';
  
  
