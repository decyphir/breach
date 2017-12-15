function params = SplotPts(P,proj,ipts,opt)
%SPLOTPTS plots parameters points (aka: the field pts) in a parameter set
% 
% Synopsis: SplotPts(P [, proj, ipts, opt])
% 
% Inputs:
%  -  P    : A Parameter set or an array of parameter set
%  -  proj : (Optional) chooses the parameters to plot; can be numbers or
%            parameters names; if not provided, is defined by
%            P(1).plot_proj if it exist. If not, is defined by the three
%            first parameter described by P(1).dim . If it contains more
%            than three elements, only the three firsts are considered.
%  -  ipts : indices of the pts to plot; all if absent or []
%  -  opt  : plot options; uses the plotting options defined in field
%            X0plot_opt or default if this is absent
% 
% see demoVDP1 for examples
% 

figure(gcf); 

if(~exist('ipts','var')||isempty(ipts))
    ipts=1:size(P(1).pts,2);
end

% dealing with plot options

if ~exist('opt','var')
    opt = [];
    if isfield(P(1),'X0plot_opt')
        opt.plot_opt = P(1).X0plot_opt;
    else
        opt.plot_opt = {'om','MarkerSize',6};
    end
elseif ischar(opt)
    opt = {opt};
end
if iscell(opt)
    opt_tmp = opt;
    opt = [];
    opt.plot_opt = opt_tmp;
end
% default plot option
if ~isfield(opt, 'plot_opt')
    opt.plot_opt = {'om','MarkerSize',6};
end

% rescaling axis

rescale = 1;
if isfield(opt, 'rescale')
    rescale = opt.rescale;
end

if(~exist('proj','var')||isempty(proj))
    if isfield(P(1),'plot_proj')
        proj = P(1).plot_proj;
    else
        switch (numel(P(1).dim))
            case {1}
                proj = P(1).dim(1);
            case {2}
                proj = P(1).dim(1: 2);
            otherwise
                proj = P(1).dim(1:3);
        end
    end
end
if ischar(proj)
    proj = {proj};
end
if ~isnumeric(proj)
    stproj = proj;
    proj = zeros(1,numel(stproj));
    for ii = 1:numel(stproj)
        proj(ii) = FindParam(P(1),stproj{ii});
    end
end
proj = proj(proj>0);
proj = proj(proj<=size(P(1).pts,1));

%  Deal with several Sampling sets

if numel(P)>1
    hold on;
    nb = numel(P);
    if nb==2
        colors = [1 0 0; 0 1 0];
    else
       colors = hsv(nb);
    end
    for ii=1:nb
        opt.plot_opt = {'+','MarkerSize',6,'Color',colors(ii,:)};
        SplotPts(P(ii),proj,[],opt);
    end
    hold off
    return
end

switch (numel(proj))
    
    case {1}
        hold on;
        params = {P.ParamList{proj(1)}};
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
        end
        ylabel('');
        set(gca, 'YTick', []);
        x = P.pts(proj(1),ipts)*rescale;
        plot(x,0*x,opt.plot_opt{:});
        
    case {2}
        hold on;
        grid on;
        params = {P.ParamList{proj(1)},P.ParamList{proj(2)}};
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
            ylabel(P.ParamList{proj(2)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
            ylabel(['x_' num2str(proj(2))],'Interpreter','tex');
        end
        
        x = P.pts(proj(1),ipts)*rescale;
        y = P.pts(proj(2),ipts)*rescale;
        plot(x,y,opt.plot_opt{:});
        
    otherwise
        hold on;
        grid on;
        params = {P.ParamList{proj(1)},P.ParamList{proj(2)},P.ParamList{proj(3)}};
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
            ylabel(P.ParamList{proj(2)},'Interpreter','none');
            zlabel(P.ParamList{proj(3)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))],'Interpreter','tex');
            ylabel(['x_' num2str(proj(2))],'Interpreter','tex');
            zlabel(['x_' num2str(proj(3))],'Interpreter','tex');
        end
        
        x = P.pts(proj(1),ipts)*rescale;
        y = P.pts(proj(2),ipts)*rescale;
        z = P.pts(proj(3),ipts)*rescale;
        plot3(x,y,z,opt.plot_opt{:});
        
end
grid on;
hold off;
%  set(gca,'FontSize',14,'FontName','times');            

end

