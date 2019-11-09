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

old_pos = get(h,'Position');
%pos([1 2]) = [0 0];
pos = old_pos;
pos(3) = old_pos(3)*x_scal;
pos(4) = old_pos(4)*y_scal;

hpos_old = old_pos(2)+old_pos(4); % pixel h position
pos(2) = hpos_old-pos(4); %makes sure top position is preserved (hopefully)
                   
set(h, 'Position', pos);
figure(h);

end
