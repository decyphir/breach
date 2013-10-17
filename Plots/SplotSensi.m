function SplotSensi(P,iX,iP,ipts)
%SPLOTSENSI Plots trajectories sensitivities
%
% Synopsis: SplotSensi(P,iX,iP,ipt)
%
% Plots trajectories sensitivities of state variables iX
% w.r.t. parameters iP
%
% Note: Uses the plotting options defined in field traj_plot_opt, and
% project on dimensions specified by field 'plot_proj'.
%
% Prerequisite: S has a traj field with an XS field. Else what's the point ?
%
% Inputs:
%   -  P        rectangular sampling set.
%
%   -  iX       indice of the X variables for which to plot the sensitivity (idem)
%
%   -  iP       indice of the sensitive parameter (optional, all if absent)
%
%   -  ipts     (optional) indices of the trajectories in S to plot (all
%               if absent)
%



%   Check inputs
if isfield(P, 'time_mult')
    time_mult = P.time_mult;
else
    time_mult = 1;
end


if isempty(P.pts)
    warning('SplotSensi:ptsEmpty','The field P.pts is empty.');
    return
end

if ~isfield(P,'traj')
    warning('SplotSensi:noTrajField','There is no field traj.');
end

if(nargin == 1)
    iX = 1:P.DimX;
    iP =P.dim;
    ipts = 1:numel(P.traj);
elseif(nargin == 2)
    iP = P.dim;
    ipts = 1:numel(P.traj);
elseif(nargin == 3)
    ipts = 1:numel(P.traj);
end

if ~isempty(iX)
    if ~isnumeric(iX)
        NiX = iX;
        iX = [];
        for ii = 1:numel(NiX)
            ind = FindParam(P,NiX{ii});
            iX(ii) = ind;
        end
    end
else
    iX = 1:P.DimX;
end
iX = iX(iX<=P.DimX);

if isempty(iP)
    iP = P.dim;
elseif ~isnumeric(iP)
    NiP = iP;
    iP = [];
    for ii = 1:numel(NiP)
        ind = FindParam(P,NiP{ii});
        iP(ii) = ind;
    end
end


if isfield(P,'traj_plot_opt')
    opt = P.traj_plot_opt;
else
    opt = {'b','LineWidth',1};
end

colors = hsv(numel(iP));

if isfield(P,'plot_proj')
    proj = P.plot_proj;
else
    proj = 1:size(P.pts,1);
end

if isfield(P,'traj_plot_opt')
    opt = P.traj_plot_opt;
end

for ii = ipts
    
    time = P.traj(ii).time;
    
    for jj = 1:numel(iX)
        if (numel(iX)>1)
            subplot(numel(iX),1,jj)
        end
        hold on;
        ylabel(['Sensi(' P.ParamList{iX(jj)} ')'],'Interpreter','none');
        
        leg = {};
        
        for k = 1:numel(iP)
            is = (find(P.dim==iP(k))-1)*size(P.traj(ii).X,1)+iX(jj);
            x = P.traj(ii).XS(is,:);
            plot(time*time_mult, x,opt{:},'Color',colors(k,:));
            leg = {leg{:}, P.ParamList{iP(k)}};
            
        end
        legend(leg{:});
    end
end

xlabel('time','Interpreter','none');

end
