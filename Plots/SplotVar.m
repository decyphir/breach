function hplot = SplotVar(P, i_var, ipts, opt, same_axe)
%SPLOTVAR plots trajectories variables separatly in the current figure if
% it exists. Otherwize, it creates a new figure.
% 
% Synopsis: hplot = SplotVar(P,[ i_var[, ipts[, opt[, same_axe]]]])
% 
% Inputs:
%  - P        : Parameter set. The trajectories must be computed or an
%               error is thrown.
%  - i_var    : (Optional, default or empty=all variables) indices or names
%               of the variables to plot
%  - ipts     : (Optional, default or empty=all trajectories) indices of
%               the parameter vectors in P for which the trajectory must be
%               plotted. For traces sets, it must set to [], so all traces
%               will be plotted (otherwise, some traces may be not
%               plotted).
%  - opt      : (Optional, default=empty) cell array containing plotting
%               options. If not defined or empty and if it exists a field
%               traj_plot_opt in P, this field is considered for opt
%  - same_axe : (Optional, default=0) boolean indicating if all variables
%               must be plotted on the same axe.
% 
% Output:
%  - hplot : array of size numel(i_var) x numel(ipts) containing handles of
%            each plotted line.
% 
% Example (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys,'x0', [-5, 5], 4);
%   P = ComputeTraj(Sys,P,0:0.01:10);
%   figure ; SplotVar(P)
%   clf ; SplotVar(P,'x0')  % plot only x0
%   clf ; SplotVar(P,[],1)  % plot only the first trajectory
%   clf ; SplotVar(P,[],1,{'--','LineWidth',3,'Color','m'}) % plot options
%   clf ; h = SplotVar(P); % plot with one color for each parameter vector
%   for ii=1:numel(h)
%       set(h(ii),'LineWidth',3); % change line width
%   end
% 
%See also SplotBoxPts SplotPts SplotTraj SplotSensi fig_resize
%

figure(gcf);

% Check inputs

if isempty(P.pts)
    error('SplotVar:emptyPtsField','The field P.pts is empty.');
end
if ~isfield(P,'traj')
    error('SplotVar:toTrajField','The parameter set has no field traj. Please compute trajectories (see ComputeTraj).');
end

if(~exist('i_var','var') || isempty(i_var))
    i_var = 1:P.DimX;
elseif(iscell(i_var) || ischar(i_var))
    i_var = FindParam(P, i_var);
end
i_var = i_var(i_var<=P.DimX);
i_var = i_var(i_var>0);

if(~exist('ipts','var')||isempty(ipts))
    ipts = 1:numel(P.traj); % manage parameter sets as trace sets
elseif isfield(P, 'traj_ref')
    ipts = unique(P.traj_ref(ipts));
end

if(~exist('opt','var')||isempty(opt))
    if isfield(P,'traj_plot_opt')
        opt = P.traj_plot_opt;
    else
        opt = [];
    end
end

if ~exist('same_axe','var')
    same_axe = 0;
end


% Plot options

if isfield(P, 'time_mult')
    time_mult = P.time_mult;
else
    time_mult = 1;
end

colors = hsv(numel(ipts));
colors = colors(:,[3 2 1]); %NM: why ?

hplot = zeros(numel(i_var),numel(ipts)); % lines handles

if(same_axe==1)
    for ii = ipts
        time = P.traj{ii}.time;
        grid on;
        %set(gca,'FontSize',12,'FontName','times');
        hold on;
        
        X = P.traj{ii}.X(i_var(:),:);
        hplot(:,ii) = plot(time*time_mult,X);
        
    end
    
    lg = cell(1,numel(i_var));
    for ii = 1:numel(i_var)
        if isfield(P,'ParamList')
            lg{ii} = P.ParamList{i_var(ii)};
        else
            lg{ii} = sprintf('x_%d',i_var(ii));
        end
    end
    hl = legend(lg);
    set(hl, 'Interpreter','none');
    hold off;
    xlabel('time')
    grid on;
else % plots on multi axes
    
    for ii = 1:numel(i_var) % for each variable
        h = subplot(numel(i_var),1,ii); % preparing the sub-graph
        grid on;
        %set(gca,'FontSize',12,'FontName','times');
        hold on;

        if isfield(P,'ParamList')
            ylabel(h,P.ParamList{i_var(ii)},'Interpreter','none');
        else
            ylabel(h, sprintf('x_%d',i_var(ii)), 'Interpreter', 'tex');
        end
        
        for jj = 1:numel(ipts)
            i_pt = ipts(jj);
            time = P.traj{i_pt}.time;
            x = P.traj{i_pt}.X(i_var(ii),:);
            if isempty(opt)
                %hplot(ii,jj) = plot(h,time*time_mult,x,'Color',colors(jj,:)); % plotting all in once with cell arrays is not faster
                hplot(ii,jj) = line(time*time_mult,x,'Color',colors(jj,:));
            else
                hplot(ii,jj) = plot(h,time*time_mult,x,opt{:});
            end
        end
        
    end
    xlabel('time')
end

end
