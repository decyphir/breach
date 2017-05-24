function [XI, YI, ZI] = QuickMesh(x,y,Z,epsi)
% 

xrange= min(x):epsi(1)/2:max(x);
yrange= min(y):epsi(2)/2:max(y);

[XI, YI] = meshgrid(xrange,yrange);
ZI = griddata(x,y,Z,XI,YI, 'linear');

surf(XI,YI,ZI,'EdgeColor','None');

