function fig_resize(f, x, y)
  
  pos = get(f,'Position');
  pos([1 2]) = [0 0];
  pos(3) = pos(3)*x;
  pos(4) = pos(4)*y;
  
  set(f, 'Position',pos);
  figure(f);