function h=flat_rect(i,d,c)

%RECT function to draw a 2-D rectangle, non transparent version
%
%Usage
%   rect(start,size,color,alpha);
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
case 0,1
    disp('Too few arguements for rect');
    return;
case 2
    l=1;    %default length of side of rectangle is 1
    c='k';  %default color of voxel is blue
 case 3
  % ok
 otherwise
    disp('Too many arguements for voxel');
end;

x=[i(1)+[0 d(1) d(1) 0]; ...
        i(2)+[0 0 d(2) d(2)]];

p = patch(x(1,:), x(2,:),'r');
set(p, 'FaceColor', 'none', 'EdgeColor',c, 'LineWidth', 2);
    

