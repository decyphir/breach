function [XI, YI, ZI] = QuickMeshSf(S,f,opt)

X = S.pts;
XI=[];
YI=[];


if (isa(f,'numeric'))
  Z = f;
else
  Z = feval(f,X);
end

switch numel(S.dim)
 case {1}
  
  if (nargin == 3)
    plot(X,Z,opt)
  else
    plot(X,Z)
  end

 case {2}
  
  x= S.pts(S.dim(1),:);
  y= S.pts(S.dim(2),:)';
  
  epsi = [min(diff(unique(sort(x)))); min(diff(unique(sort(y))))];
  [XI, YI, ZI] = QuickMesh(x,y,Z,epsi);
  
 otherwise
  disp(['we should now enter the forth dimension, tin din tin din' ...
	' tin din '])
  
end
