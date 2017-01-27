function fig_resize(h, x_scal, y_scal)
%FIG_RESIZE is a simple function to resize a Matlab figure by a vertical
% and horizont scaling factor
% 
% Synopsis: fig_resize(h, x_scal, y_scal)
% 
% Inputs:
%  - h      : handle of the figure to rescale
%  - x_scal : scale on x
%  - y_scal : scale on y
% 
% Example (Lorentz84):
%   CreateSystem
%   P = CreateParamSet(Sys,'x0');
%   P = ComputeTraj(Sys,P,0:0.1:10);
%   SplotVar(P)
%   fig_resize(gcf, 1, 2) % multiplies the vertical size by 2
%   fprintf(' press Enter to continue');
%   pause
%   fig_resize(gcf, 2, 0.5)
%

pos = get(h,'Position');
pos([1 2]) = [0 0];
pos(3) = pos(3)*x_scal;
pos(4) = pos(4)*y_scal;

set(h, 'Position', pos);
figure(h);

end
