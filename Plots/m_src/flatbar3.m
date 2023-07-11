function flatbar3(X,Y,Z,widthx,widthy)

[n_y, n_x]=size(Z);

if isscalar(widthx)
    widthx = ones(1,n_x)*widthx;
end

if isscalar(widthy)
    widthy = ones(1,n_y)*widthy;
end

for j=1:n_y
    for k=1:n_x
        if ~isnan(Z(j,k))
            drawbar(X(j,k),Y(j,k),Z(j,k),widthx(k)/2,widthy(j)/2)
        end
    end
end

zlim=[min(Z(:)) max(Z(:))];
if zlim(1)>0
    zlim(1)=0;
end
if zlim(2)<0
    zlim(2)=0;
end

axis([min(X(:))-widthx(1) max(X(:))+widthx(end) min(Y(:))-widthy(1) max(Y(:))+widthy(end) zlim])
Zmin = min(Z(:));
Zmax = max(Z(:));
if Zmin==Zmax
    Zmax = Zmin+min(abs(Zmax),1)/1e3;
end
caxis([Zmin Zmax]);

end

function drawbar(x,y,z,widthx,widthy)
h(1)=patch([-widthx -widthx widthx widthx]+x,[-widthy widthy widthy -widthy]+y,[0 0 0 0],'b');
set(h,'facecolor','flat','FaceVertexCData',z)
end