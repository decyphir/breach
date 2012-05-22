function QuickContour(x,y,Z,epsi, val)
 
  if ~exist('val')
    val = [0. 0.];
  end
  xrange= min(x)-epsi(1)/2:epsi(1)/2:max(x)+epsi(1)/2;
  yrange= min(y)-epsi(2)/2:epsi(2)/2:max(y)+epsi(2)/2;

  [XI YI] = meshgrid(xrange,yrange);
  ZI = griddata(x,y,Z,XI,YI,[], {'Qt', 'Qbb','Qc', 'QbB'});

  contour(XI,YI,ZI,val);
