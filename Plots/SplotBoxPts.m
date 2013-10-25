function SplotBoxPts(P, proj, ipts, opt, col, alph)
%SPLOTBOXPTS  Plots the points in the field pts of a parameter set and the
% boxes they represent.
%
% Synopsis:  SplotBoxPts(P, proj, ipts, opt, col, alpha)
%
% Inputs:
%  -  P        parameter  set.
%  -  proj     chooses the parameters to plot; can be numbers or
%              parameters names; all if []
%  -  ipts     indices of the pts to plot; all if absent or []
%  -  opt      Uses the plotting options defined in field X0plot_opt or
%              default if this is absent
%  -  col      Color (e.g. 'g' 'r' 'k' etc)
%  -  alph     alpha rendering
% 
% Outputs:
%  - none, but a figure
%

if ~exist('proj','var')
    proj = [];
end
if(isempty(proj) && isfield(P,'plot_proj'))
    proj = P.plot_proj;
end
if ~isnumeric(proj)
    proj = FindParam(P,proj);
end
proj = proj(proj>0);
proj = proj(proj<=size(P.pts,1));
if isempty(proj)
    proj = 1:min(3,numel(P.dim));
    proj = P.dim(proj);
elseif(numel(proj)>3)
    proj = proj(1:3);
end

if(~exist('ipts','var')||isempty(ipts))
    ipts = 1:size(P.pts,2);
end

if(exist('opt','var')&&(~isempty(opt)))
    SplotPts(P, proj, ipts, opt);
else
    SplotPts(P, proj, ipts);
end

nb_pts = numel(P.pts(1,:));
if(nb_pts==0)
    return;
end

if ~exist('col','var')
    col = 'b';
end

if ~exist('alph','var')
    alph = .1;
end


hold on;

switch(numel(proj))
    
    case 1
        
        hold on;
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
        end
        
        DX(2) = 0;
        i = find(P.dim == proj);
        
        if ~isempty(i)
            
            for k=ipts
                
                DX(1) = P.epsi(i,k);
                X = [P.pts(proj,k)'-DX(1),-DX(2)];
                rect(X,2*DX,col,alph);
                
            end
        end
        %set(gca, 'YLim', [-1 1], 'YtickLabel', {},'Interpreter','none');
        
    case 2
        nb_pts = size(P.pts,2);
        
        for k=ipts
            for j = 1:2
                i = find(P.dim == proj(j));
                if isempty(i)
                    DX(j) = 0;
                else
                    DX(j) = P.epsi(i,k);
                end
            end
            X = P.pts(proj,k)'-DX;
            rect(X,2*DX,col,alph);
        end
        
    case 3
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
            ylabel(P.ParamList{proj(2)},'Interpreter','none');
            zlabel(P.ParamList{proj(3)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
            ylabel(['x_' num2str(proj(2))],'Interpreter','tex');
            zlabel(['x_' num2str(proj(3))],'Interpreter','tex');
        end
        nb_pts = size(P.pts,2);
        
        for k=ipts
            for j = 1:3
                i = find(P.dim == proj(j));
                if (isempty(i)||P.epsi(i,k)==0)
                    DX(j) = 0;
                else
                    DX(j) = P.epsi(i,k);
                end
            end
            
            X = P.pts(proj,k)'-DX;
            voxel(X,2*DX,col,alph);
        end
end
grid on;
hold off;

end
