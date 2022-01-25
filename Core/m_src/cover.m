function [cnt, grd] = cover(x, delta)
% Compute grid of interest 
xmin = min(x); 
xmax = max(x);

if nargin==1  % sth like 10 bins by default..
    delta = (xmax-xmin)/10; 
end

grd = xmin:delta:xmax;

Glow = repmat(grd', 1, size(x,2));
Gup = Glow(2:end,:);
Gup(end+1,:) = inf;

X = repmat(x, size(grd,2),1);
Xin = (X==Glow)|(X>Glow)&(X<Gup)|X==Gup;

cnt = sum(Xin,2)';


