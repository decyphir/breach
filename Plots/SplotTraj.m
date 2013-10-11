function SplotTraj(P,proj,ipts,opt,t0)
%SPLOTTRAJ plots the trajectories in field traj of P in the phase space.
%
%   Note:  Uses the plotting options defined in field traj_plot_opt, and
%   project on dimensions specified by field 'plot_proj' if these fields
%   are defined
%
% Synopsis: SplotTraj(S [, proj, ipts, opt, t0])
%
%   Inputs:
%    -  P    : Parameter set. It must have a traj field, else what's the
%              point (and an error is thrown)?
%    -  proj : (optional) variables to plot (the three first if [])
%    -  ipts : indices of the initial pts from which to plot
%    -  opt  : plotting option e.g: {'r','LineWidth',4}
%    -  t0   : (optional) starting time to plot trajectories from, zero by
%              default.
%
%   Outputs:
%    - None, but plot a figure, hopefully an interesting one.
%

% Check inputs
if isempty(P.pts)
    error('SPlotTraj:noPts','P.pts is empty !');
end

if ~isfield(P,'traj')
    error('SPlotTraj:noTraj','No trajectory computed for this set')
end

if isfield(P, 'time_mult')
    time_mult = P.time_mult;
else
    if isfield(P,'opt')
        if isfield(P.opt, 'time_mult')
            time_mult= P.opt.time_mult;
        else
            time_mult=1;
        end
    else
        time_mult=1;
    end
end

if isfield(P, 'rescale')
    rescale = P.rescale;
else
    if isfield(P,'opt')
        if isfield(P.opt, 'rescale')
            rescale = P.opt.rescale;
        else
            rescale = 1;
        end
    else
        rescale = 1;
    end
end


if ~exist('t0','var')
    t0 = 0;
end


if isfield(P,'traj_plot_opt')
    opt = P.traj_plot_opt;
elseif(~exist('opt','var')||isempty(opt))
    opt = {'b','LineWidth',1};
end

if ischar(opt)
    opt = {opt};
end


% find the projected axes
if isfield(P,'plot_proj')
    proj = P.plot_proj;
elseif(~exist('proj','var')||isempty(proj))
    switch (P.DimX)
        case {1}
            proj=1;
        case {2}
            proj=[ 1 2 ];
        otherwise
            proj = [ 1 2 3 ];
    end
end
if ischar(proj)
    proj = {proj};
end
if ~isnumeric(proj)
    stproj = proj;
    proj=zeros(1,numel(stproj));
    for ii = 1:numel(stproj)
        proj(ii) = FindParam(P,stproj{ii});
    end
end
proj = proj(proj~=0);


if(~exist('ipts','var')||isempty(ipts))
    ipts = 1:numel(P.traj);
end
if isfield(P, 'traj_ref')
    ipts = unique(P.traj_ref(ipts));
end
ipts = ipts(ipts>0); % to avoid non computed trajectories

switch (numel(proj))
    
    case {1}
        hold on;
        xlabel('time')
        
        if isfield(P,'ParamList')
            ylabel(P.ParamList{proj(1)},'Interpreter','none');
        else
            ylabel(['x_' num2str(proj(1))]);
        end
        
        if isfield(P,'traj_plot_opt')
            opt = P.traj_plot_opt;
        end
        
        for ii = ipts
            time = P.traj(ii).time;
            if(numel(time>t0))
                y = P.traj(ii).X(proj(1),time>=t0) * rescale;
                x = P.traj(ii).time(time>=t0)*time_mult;
                plot(x,y, opt{:});
                % drawnow
            end
        end
        %drawnow
        
    case {2}
        hold on;
        
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
            ylabel(P.ParamList{proj(2)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))]);
            ylabel(['x_' num2str(proj(2))]);
        end
        
        if isfield(P,'traj_plot_opt')
            opt = P.traj_plot_opt;
        end
        
        for ii = ipts
            time = P.traj(ii).time;
            if (numel(P.traj(ii).time)>0)
                x = P.traj(ii).X(proj(1),time>=t0)*rescale;
                y = P.traj(ii).X(proj(2),time>=t0)*rescale;
                plot(x,y,opt{:});
            end
        end
        %drawnow
        
    otherwise
        hold on;
        
        if isfield(P,'ParamList')
            xlabel(P.ParamList{proj(1)},'Interpreter','none');
            ylabel(P.ParamList{proj(2)},'Interpreter','none');
            zlabel(P.ParamList{proj(3)},'Interpreter','none');
        else
            xlabel(['x_' num2str(proj(1))]);
            ylabel(['x_' num2str(proj(2))]);
            zlabel(['x_' num2str(proj(3))]);
        end
        
        if isfield(P,'traj_plot_opt')
            opt = P.traj_plot_opt;
        end
        
        for ii = ipts
            time = P.traj(ii).time;
            if ~isempty(time)
                x = P.traj(ii).X(proj(1),time>=t0)*rescale;
                y = P.traj(ii).X(proj(2),time>=t0)*rescale;
                z = P.traj(ii).X(proj(3),time>=t0)*rescale;
                plot3(x,y,z,opt{:});
            end
            %drawnow
        end
end

grid on;
if(t0==0)&&(numel(proj)>1)
    SplotPts(P,proj,ipts);
    if isfield(P,'Xf')
        SplotXf(P,proj,ipts);
    end
end

end
