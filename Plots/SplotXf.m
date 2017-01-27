function SplotXf(P, proj, ipts, opt)
%SPLOTXF plots in the current figure final points of computed trajectories
% (field Xf) of a parameter set.
% 
% Synopsis: SplotXf(P [, proj, ipts, opt])
% 
% Inputs:
%  -  P    : Parameter set. Any set with a field 'Xf' will do though.
%  -  proj : (Optional, default or empty=P.plot_proj if exists, three first
%            variable of P otherwise) chooses the parameters to plot. It
%            can be numbers or parameters names.
%  -  ipts : (Optional, default or empty=all parameter vectors) choses the
%            indexes of the parameter vectors to plot.
%  -  opt  : Uses the plotting options defined in field Xfplot_opt or
%            default if this is absent
% 
% Output:
%  - none, but a figure
% 
% Example (Lorentz84):
%  CreateSystem;
%  P = CreateParamSet(Sys,{'F','G'},[10,20;0,4],4); % Create 16 parameter vectors
%  P = ComputeTraj(Sys,P,{0,10});
%  figure ; SplotXf(P);
%  SplotTraj(P); % shows lines from (0,0,0), which is the origine
% 
%See also ComputeTraj
%


if(~exist('ipts','var')||isempty(ipts))
    ipts = 1:size(P.pts,2);
end
if isfield(P,'traj_ref')
    ipts = unique(P.traj_ref(ipts));
end
ipts = ipts(ipts>0); % to avoid non-computed traj

if ~exist('opt','var')
    opt = [];
    if isfield(P,'X0plot_opt')
        opt.plot_opt = P.X0plot_opt;
    else
        opt.plot_opt = {'+k','MarkerSize',6};
    end
else
    if ischar(opt)
        opt = {opt};
    end
end

if iscell(opt)
    opt0 = opt;
    opt = [];
    opt.plot_opt = opt0;
end

% default plot option
if ~isfield(opt, 'plot_opt')
    opt.plot_opt = {'+k','MarkerSize',6};
end

% rescaling axis

rescale = 1;

if isfield(opt, 'rescale')
    rescale = opt.rescale;
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
proj = proj(proj<=P.DimX);
if isempty(proj)
    proj = 1:min(3,P.DimX);
elseif(numel(proj)>3)
    proj = proj(1:3);
end



switch(numel(proj))
    case {1}
        hold on;
        x = P.Xf(proj(1),ipts);
        plot(x,0*x,opt.plot_opt{:});
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))]);
        end
        
    case {2}
        hold on;
        x = P.Xf(proj(1),ipts);
        y = P.Xf(proj(2),ipts);
        plot(x,y,opt.plot_opt{:});
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
            ylabel(P.ParamList{proj(2)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))]);
            ylabel(['x_' num2str(proj(2))]);
        end
        
    otherwise
        hold on;
        x = P.Xf(proj(1),ipts);
        y = P.Xf(proj(2),ipts);
        z = P.Xf(proj(3),ipts);
        plot3(x,y,z,opt.plot_opt{:});
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
            ylabel(P.ParamList{proj(2)},'Interpreter','none');
            zlabel(P.ParamList{proj(3)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))]);
            ylabel(['x_' num2str(proj(2))]);
            zlabel(['x_' num2str(proj(3))]);
        end
        
end

grid on;

end
