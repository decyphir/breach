function Sh = SobolRefine(S, nb)
% QUASIREFINE  Sample quasi-uniformly a parameter set 
% 
% Synopsis:  Ph = QuasiRefine(S, nb)
%  
% Example: 
%
%   CreateSystem;
%   P = CreateSampling(Sys); % Create default parameter set for system Sys
%   Ph = SobolRefine(P, 1000); % Sample with 1000 points
%  
%   SplotBoxPts(P); % Parameter set before sampling 
%   SplotPts(Ph);   % plots the generated points
%
% TODO optionalize the sequence used - so far Sobol
%  
  
  Sh = SobolRefine(S,nb);