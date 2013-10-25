function SplotSensi(P, iX, iP, ipts)
%SPLOTSENSI Plots trajectories sensitivities of state variables iX w.r.t.
% parameters iP.
%
% Synopsis: SplotSensi(P, iX, iP, ipts)
%
% Note: Uses the plotting options defined in field traj_plot_opt, and
% project on dimensions specified by field 'plot_proj'.
%
% Inputs:
%  -  P    : Parameter set. Trajectories and the sensitivities must be
%            computed. If the trajectories are not computed, an error is
%            thrown.
%  -  iX   : (Optional, default=all variables) indice or names of the
%            variables for which the sensitivity is plotted
%  -  iP   : (Optional, default=P.dim) indexes or names of the sensitive
%            parameter
%  -  ipts : (Optional, default=all traj) indices of the trajectories in P
%            to plot
%


if isempty(P.pts)
    error('SplotSensi:ptsEmpty','The field P.pts is empty.');
end
if ~isfield(P,'traj')
    error('SplotSensi:noTrajField','There is no field traj.');
end
if ~isfield(P,'type')
    P.type = '';
end
if(strcmp(P.type,'') && ~isfield(P,'traj_ref'))
    P.traj_ref = 1:size(P.pts,2);
end

if ~exist('iX','var')
    iX = [];
elseif ~isnumeric(iX)
    iX = FindParam(iX);
end
iX = iX(iX>0);
iX = iX(iX<=P.DimX);
if isempty(iX)
    iX = 1:P.DimX;
end

if ~exist('iP','var')
    iP = [];
elseif ~isnumeric(iP)
    iP = FindParam(iP);
end
iP = iP(iP>0);
iP = iP(iP<=size(P.pts,1));
if isempty(iP)
    iP = P.dim;
end

if ~exist('ipts','var')
    ipts = [];
end
if strcmp(P.type,'')
    ipts = ipts(ipts>0);
    ipts = ipts(ipts<=size(P.pts,2));
    if isempty(ipts)
        ipts = unique(P.traj_ref(1:size(P.pts,2)));
    end
    ipts = ipts(ipts~=0); % avoid not computed trajectories
else
    ipts = ipts(ipts>0);
    ipts = ipts(ipts<=numel(P.traj));
    if isempty(ipts)
        ipts = 1:numel(P.traj);
    end
end


%   Check inputs
if isfield(P, 'time_mult')
    time_mult = P.time_mult;
else
    time_mult = 1;
end

if isfield(P,'traj_plot_opt')
    opt = P.traj_plot_opt;
else
    opt = {'b','LineWidth',1};
end

colors = hsv(numel(iP));

if isfield(P,'traj_plot_opt')
    opt = P.traj_plot_opt;
end

for jj = 1:numel(iX) % we prepare the graphic
    subplot(numel(iX),1,jj);
    hold on;
    grid on;
    ylabel(['Sensi(' P.ParamList{iX(jj)} ')'],'Interpreter','none');
end

for ii = ipts % then, we plot
    time = P.traj(ii).time;
    
    for jj = 1:numel(iX)
        subplot(numel(iX),1,jj)
        for kk = 1:numel(iP)
            is = (find(P.dim==iP(kk))-1)*size(P.traj(ii).X,1)+iX(jj);
            x = P.traj(ii).XS(is,:);
            plot(time*time_mult, x, opt{:}, 'Color', colors(kk,:));
        end
    end
end

leg = P.ParamList(iP); % we show the legend.
for jj = 1:numel(iX)
    subplot(numel(iX),1,jj);
    legend(leg);            %%%%% THIS LINE TAKES A VERY LONG TIME !!! %%%%
end
xlabel('time','Interpreter','none');

end
