% function f = cam(x)
% Four-hump camel
function f = cam(x)
f=(4-2.1.*x(1).^2+x(1).^4./3).*x(1).^2+x(1).*x(2)+(-4+4.*x(2).^2).*x(2).^2;        

