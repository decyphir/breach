function h=rect(i,d,c,alpha)
%RECT function to draw a 2-D rectangle
%
% Synopsis: rect(start, size, color, alpha);
%
%   will draw a rectangle at 'start' of size 'size' of color 'color' and
%   transparency alpha (1 for opaque, 0 for transparent)
%   Default size is 1
%   Default color is blue
%   Default alpha value is 1
%
%   start is a three element vector [x,y]
%   size the a three element vector [dx,dy]
%   color is a character string to specify color
%       (type 'help plot' to see list of valid colors)
%
%
%   rect([2 3] ,[1 2],'r',0.7);
%   axis([0 10 0 10]);


switch(nargin),
    case 0
        disp('Too few arguements for rect');
        return;
    case 1
        l=1;    %default length of side of rectangle is 1
        c='b';  %default color of voxel is blue
        alpha =1;
    case 2,
        c='b';
        alpha=1;
    case 3,
        alpha=1;
    case 4,
        %do nothing
    otherwise
        disp('Too many arguements for voxel');
end;

x=[i(1)+[0 d(1) d(1) 0]; ...
    i(2)+[0 0 d(2) d(2)]];
if strcmp(c,'none')
    p = patch(x(1,:), x(2,:),'r');
    set(p, 'FaceColor', 'none');
    return
end
h = fill(x(1,:), x(2,:),c,'FaceAlpha',alpha);

end
