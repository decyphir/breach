function [cnt, grd] = cover(x, delta)
% Compute grid of interest 
xmin = min(x); 
xmax = max(x);

if nargin==1  % sth like 10 bins by default..
    delta = (xmax-xmin)/10; 
end

grd = xmin:delta:xmax;

Gabove = repmat(grd', 1, size(x,2));
Gbelow = Gabove(2:end,:);
Gbelow(end+1,:) = inf;

X = repmat(x, size(grd,2),1);
Xin = (X>= Gabove)&(X<Gbelow);

cnt = sum(Xin,2)';
