function QuickMesh(x,y,Z,epsi)

xrange= min(x)-epsi(1)/2:epsi(1)/2:max(x)+epsi(1)/2;
yrange= min(y)-epsi(2)/2:epsi(2)/2:max(y)+epsi(2)/2;

[XI YI] = meshgrid(xrange,yrange);
ZI = griddata(x,y,Z,XI,YI, 'linear', {'Qt', 'Qbb','Qc', 'QbB'} );
%ZI = griddata(x,y,Z,XI,YI, 'cubic', {'Qt', 'Qbb','Qc', 'QbB'} );
%mesh(XI,YI,ZI);

surf(XI,YI,ZI,'EdgeColor','None');
%contour(XI,YI,ZI,[.0 .0]);
