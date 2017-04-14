function SplotBoxPts(P, proj, ipts, opt, col, alph)
%SPLOTBOXPTS plots the points in the field pts of a parameter set and the
% boxes they represent.
% 
% Synopsis:  SplotBoxPts(P, proj, ipts, opt, col, alpha)
% 
% Inputs:
%  -  P    : parameter  set. It may contain many parameter vectors
%  -  proj : (Optional, default of empty=P.plot_proj if exist or
%            P.dim[1:3]) chooses the parameters to plot. It can be indexes
%            or parameters names
%  -  ipts : (Optional, default or empty=all parameter vectors) indices of
%            the parameter vectors to plot
%  -  opt  : (Optional, default=empty) uses the plotting options defined in
%            field X0plot_opt or default if this is absent
%  -  col  : (Optional, default='b') caracter indicating the color for the
%            plot (e.g. 'g' 'r' 'k' etc)
%  -  alph : (Optional, default=0.1) alpha rendering
% 
% Output:
%  - none, but a figure
% 
% Example (Lorentz84):
%   See SConcat
%

if(~isfield(P,'pts') || size(P.pts,2)==0)
    return;
end

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

if ~exist('col','var')
    col = 'b';
end

if ~exist('alph','var')
    alph = .1;
end

hold on;

switch(numel(proj))
    
    case 1
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
        end
        
        DX(2) = 0;
        idx = find(P.dim == proj);
        
        for kk=ipts % for each parameter vector
            if ~isempty(idx)
                DX(1) = P.epsi(idx,kk);
            else
                DX(1) = 0;
            end
            X = [P.pts(proj,kk)-DX(1),-DX(2)];
            rect(X, 2*DX, col, alph);
        end
        %set(gca, 'YLim', [-1 1], 'YtickLabel', {},'Interpreter','none');
        
    case 2
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
            ylabel(P.ParamList{proj(2)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
            ylabel(['x_' num2str(proj(2))],'Interpreter','tex');
        end
        
        DX = zeros(1,2);
        for kk=ipts
            for jj = 1:2 % if proj(jj) in P.dim, set the width
                idx = find(P.dim == proj(jj));
                if isempty(idx)
                    DX(jj) = 0;
                else
                    DX(jj) = P.epsi(idx,kk);
                end
            end
            X = P.pts(proj,kk)'-DX;
            rect(X, 2*DX, col, alph);
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
        
        DX = zeros(1,3);
        for kk=ipts
            for jj = 1:3 % if proj(jj) in P.dim, set the width
                idx = find(P.dim == proj(jj));
                if isempty(idx)
                    DX(jj) = 0;
                else
                    DX(jj) = P.epsi(idx,kk);
                end
            end
            
            X = P.pts(proj,kk)'-DX;
            voxel(X, 2*DX, col, alph);
        end
end
grid on;
hold off;

end
