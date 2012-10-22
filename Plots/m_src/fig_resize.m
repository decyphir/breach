function fig_resize(f, x, y)
% FIG_RESIZE Simple function to resize a Matlab figure by a vertical and
% horizont scaling factor
%
% Synopsys: fig_resize(handle, x_scal, y_scal)  
%
% Example: fig_resize(gcf, 1, 2) multiplies the vertical size by 2
  
  
  pos = get(f,'Position');
  pos([1 2]) = [0 0];
  pos(3) = pos(3)*x;
  pos(4) = pos(4)*y;
  
  set(f, 'Position',pos);
  figure(f);
  