function SplotVarExp(P, iX, ipts)
%SPLOTVAREXP plots trajectories states variables and represent expansion
% given by sensitivity
% 
% Synopsis: SplotVarExp(P[, iX[, ipts]])
% 
% Note:   Uses the plotting options defined in field traj_plot_opt, and
% project on dimensions specified by field 'plot_proj'.
% 
% Inputs:
%  -  P    : Parameter set. It must have a traj field.
%  -  iX   : (Optional, default or empty=all variables) indices or names of
%            the variables to plot.
%  -  ipts : (Optional, default or empty=all trajectories) indices of the
%            parameter vectors for which the trajectories are plotted.
% 
% Output:
%  -  Some figure, hopefully an interesting one.
% 
% Example (Lorentz84):
%   CreateSystem
%   P = CreateParamSet(Sys,'a',[0,0.5]);
%   P = ComputeTrajSensi(Sys,P,0:0.1:10);
%   figure ; SplotVarExp(P)  % nice picture uh!
%


% Check inputs
if isempty(P.pts)
    error('SplotVarExp:ptsEmpty','The field P.pts is empty.');
end

if ~isfield(P,'traj')
    error('SplotVarExp:noTrajField','No trajectory computed for this set.')
end

if(~exist('iX','var') || isempty(iX))
    iX = 1:P.DimX;
elseif ~isnumeric(iX)
    iX = FindParam(iX);
end
iX = iX(iX>0);
iX = iX(iX<=P.DimX);

if(~exist('ipts','var') || isempty(ipts))
    ipts = 1:numel(P.traj);
else
    ipts = unique(P.traj_ref(ipts));
end

% manage options
if isfield(P,'traj_plot_opt')
    opt = P.traj_plot_opt;
else
    opt = {'b','LineWidth',1};
end

% plot
for ii = ipts
    time = P.traj{ii}.time;
    
    for jj = 1:numel(iX)
        subplot(numel(iX),1,jj)
        hold on;

        if isfield(P,'ParamList')
            ylabel(P.ParamList{iX(jj)},'Interpreter','none');
        else
            ylabel(['x_' num2str(iX(jj))]);
        end

        x = P.traj{ii}.X(iX(jj),:);
        e = P.traj{ii}.Expa(iX(jj),:);
        base = min(x-e-.1);
        
        area(time,x+e,base,'FaceColor',[.5 .9 .6],...
            'EdgeColor','w',...
            'LineWidth',1)
        hold on
        area(time,x-e,base,'FaceColor',[1 1 1],...
            'EdgeColor','w',...
            'LineWidth',1)
        plot(time,x+e,'k');
        plot(time,x-e,'k');
        plot(time,x,opt{:});
    end
end
xlabel('time')

end
