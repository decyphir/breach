function SplotSensi(P, iX, iP, ipts)
%SPLOTSENSI plots trajectories sensitivities of state variables iX w.r.t.
% parameters iP.
% 
% Synopsis: SplotSensi(P, iX, iP, ipts)
% 
% Note: Uses the plotting options defined in field traj_plot_opt, and
% project on dimensions specified by field 'plot_proj'.
% 
% Inputs:
%  - P    : Parameter set. Trajectories and the sensitivities must be
%           computed. If the trajectories are not computed, an error is
%           thrown.
%  - iX   : (Optional, default=all variables) indice or names of the
%           variables for which the sensitivity is plotted
%  - iP   : (Optional, default=P.dim) indexes or names of the sensitive
%           parameter
%  - ipts : (Optional, default=all traj) indices of the trajectories in P
%           to plot
% 
% Output:
%  - none, but a figure.
% 
% Example (Lorentz84):
%   CreateSystem
%   P = CreateParamSet(Sys, {'x0','x1','a'}, [0,1 ; 0,1 ; 0,0.5]);
%   P = ComputeTrajSensi(Sys,P,0:0.001:0.1);
%   figure ; SplotSensi(P)
% 
%See also ComputeTrajSensi SplotVar SplotBoxPts SplotPts SplotTraj
%

% Check inputs
if isempty(P.pts)
    error('SplotSensi:ptsEmpty','The field P.pts is empty.');
end
if ~isfield(P,'traj')
    error('SplotSensi:noTrajField','There is no field traj.');
end
if ~isfield(P,'type')
    P.type = 'Breach';
end
if(strcmp(P.type,'') && ~isfield(P,'traj_ref'))
    P.traj_ref = 1:size(P.pts,2);
end

if ~exist('iX','var')
    iX = [];
elseif ~isnumeric(iX)
    iX = FindParam(P,iX);
end
iX = iX(iX>0);
iX = iX(iX<=P.DimX);
if isempty(iX)
    iX = 1:P.DimX;
end

if ~exist('iP','var')
    iP = [];
elseif ~isnumeric(iP)
    iP = FindParam(P,iP);
end
iP = iP(iP>0);
iP = iP(iP<=size(P.pts,1));
if isempty(iP)
    iP = P.dim;
end

if ~exist('ipts','var')
    ipts = [];
end
if strcmp(P.type,'Breach')
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

% manage options
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

% the plot itself
for jj = 1:numel(iX) % we prepare the graphic
    subplot(numel(iX),1,jj);
    hold on;
    grid on;
    ylabel(['Sensi(' P.ParamList{iX(jj)} ')'],'Interpreter','none');
end

for ii = 1:numel(iX) % then, we plot
    subplot(numel(iX),1,ii)
    for jj = 1:numel(ipts)
        time = P.traj{ipts(jj}).time;
        for kk = 1:numel(iP)
            is = (find(P.dim==iP(kk))-1)*size(P.traj{ipts(jj}).X,1)+iX(ii);
            x = P.traj{ipts(jj}).XS(is,:);
            plot(time*time_mult, x, opt{:}, 'Color', colors(kk,:));
        end
        if(jj==1) % we plot the legend when only the sensitivity of first parameter vector is
            legend(P.ParamList(iP)); % plot to avoid to loose time when showing the legend
        end
    end
end
xlabel('time','Interpreter','none');

end
