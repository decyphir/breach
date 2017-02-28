function [cnt, grd1,grd2] =cover2d(x, delta1, delta2)

%% x1
x1 = x(1,:);
% Compute grid of interest for dim 1
xmin1 = min(x1); 
xmax1 = max(x1);

if nargin<=1
    delta1 = (xmax1-xmin1)/10; 
end
grd1 = xmin1:delta1:xmax1;

Gabove1 = repmat(grd1', 1, size(x1,2));
Gbelow1 = Gabove1(2:end,:);
Gbelow1(end+1,:) = inf;

X1 = repmat(x1, size(grd1,2),1);
Xin1 = (X1>= Gabove1)&(X1<Gbelow1);

%% x2
x2 = x(2,:);
% Compute grid of interest for dim 2
xmin2 = min(x2); 
xmax2 = max(x2);

if nargin<=2
    delta2 = (xmax2-xmin2)/10; 
end

grd2 = xmin2:delta2:xmax2;

Gabove2 = repmat(grd2', 1, size(x2,2));
Gbelow2 = Gabove2(2:end,:);
Gbelow2(end+1,:) = inf;

X2 = repmat(x2, size(grd2,2),1);
Xin2 = (X2>= Gabove2)&(X2<Gbelow2);

%% TODO find something smarter

for i = 1:numel(grd1)
    for j= 1:numel(grd2)
       cnt(i,j) = sum( Xin1(i,:)&Xin2(j,:));      
    end
end