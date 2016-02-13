function [S  Z] = QuickMeshSf(S,f,opt)

X = S.pts;

if (isa(f,'numeric'))
    Z=f;
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
        epsi = [min(S.epsi(1,:)) ; min(S.epsi(2,:))];
        QuickContour(x,y,Z,epsi);
        
    otherwise
        disp(['we should now enter the forth dimension, tin din tin din' ...
            ' tin din '])
        
end
