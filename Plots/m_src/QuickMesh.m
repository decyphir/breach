function [XI YI ZI] = QuickMesh(x,y,Z,epsi)

%xrange= min(x)-epsi(1)/2:epsi(1)/2:max(x)+epsi(1)/2;
%yrange= min(y)-epsi(2)/2:epsi(2)/2:max(y)+epsi(2)/2;

xrange= min(x):epsi(1)/2:max(x);
yrange= min(y):epsi(2)/2:max(y);

[XI YI] = meshgrid(xrange,yrange);

ZI = griddata(x,y,Z,XI,YI, 'linear');

ex = log10(min(xrange(1)))+1;
ey = log10(min(yrange(1)))+1;
e = max(ex,ey);

if (e>=-6)
    e=1
end
surf(XI*10^(-e),YI*10^(-e),ZI,'EdgeColor','None');

