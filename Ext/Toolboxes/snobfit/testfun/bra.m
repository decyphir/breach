% function f = bra(x)
% Branin's function
function f = bra(x)
a=1;
b=5.1/(4*pi*pi);
c=5/pi;
d=6;
h=10;
ff=1/(8*pi);
f=a.*(x(2)-b.*x(1).^2+c.*x(1)-d).^2+h.*(1-ff).*cos(x(1))+h;
