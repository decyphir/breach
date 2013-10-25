function SplotVar(P, i_var, ipts, opt, same_axe)
%SPLOTVAR Plots trajectories variables separatly
%
% Synopsis:  SplotVar(P, [i_var[, ipts[, opt[, same_axe]]]])
%
% Inputs:
%  -  P       : Parameter set. The trajectories must be computed or an
%               error is thrown.
%  -  i_var   : (Optional, default or empty=all variables) indices or names
%               of the variables to plot
%  -  ipts    : (Optional, default or empty=all trajectories) indices of
%               the parameter vectors in P for which the trajectory must be
%               plotted. For traces sets, it must set to [], so all traces
%               will be plotted (otherwise, some traces may be not
%               plotted).
%  - opt      : (Optional) plotting options
%  - same_axe : (Optional, default=0) boolean indicating if all variables
%               must be plotted on the same axe.
%
% Output:
%  - none, but a figure
% 
% Example (Lorentz84):
%  CreateSystem;
%  P = CreateParamSet(Sys,'x0', [-5, 5], 4);
%  P = ComputeTraj(Sys,P,0:0.01:10);
%  figure ; SplotVar(P)
%  clf ; SplotVar(P,'x0')  % plot only x0
%  clf ; SplotVar(P,[],1)  % plot only the first trajectory
%

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
colors = colors(:,[3 2 1]);

if(same_axe==1)
    for ii = ipts
        time = P.traj(ii).time;
        grid on;
        %set(gca,'FontSize',12,'FontName','times');
        hold on;
        
        X = P.traj(ii).X(i_var(:),:);
        plot(time*time_mult,X);
        
    end
    
    lg = cell(1,numel(i_var));
    for ii = 1:numel(i_var)
        if isfield(P,'ParamList')
            lg{ii} = P.ParamList{i_var(ii)};
        else
            lg{ii} = ['x_' num2str(i_var(ii))];
        end
    end
    hl = legend(lg);
    set(hl, 'Interpreter','none');
    hold off;
    xlabel('time')
    
else % plots on multi axes
    for jj = 1:numel(i_var) % preparing the graph
        subplot(numel(i_var),1,jj)
        grid on;
        %set(gca,'FontSize',12,'FontName','times');
        hold on;

        if isfield(P,'ParamList')
            ylabel(P.ParamList{i_var(jj)},'Interpreter','none');
        else
            ylabel(['x_' num2str(i_var(jj))],'Interpreter','tex');
        end
    end
    
    ci = 1;
    for ii = ipts % then plotting
        time = P.traj(ii).time;
        
        for jj = 1:numel(i_var)
            subplot(numel(i_var),1,jj)
            
            x = P.traj(ii).X(i_var(jj),:);
            if isempty(opt)
                plot(time*time_mult, x, 'Color', colors(ci,:));
            else
                plot(time*time_mult, x, opt{:});
            end
        end
        ci = ci+1;
    end
    hold off;
    xlabel('time')
end

end
