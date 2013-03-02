function [Mend, opts] =  SplotSensiBar(Sys, P, ipts, opts)
%SPLOTSENSIBAR plots an histogram of sensitivities
%
% Synopsis:  [M, opts] = SplotSensiBar(Sys, P, ipts [, opts])
%
% Input:
%   - Sys  : the system
%   - P    : a parameter set containing a field pts
%   - ipts : indexes of the parameter sets of P 
%   - opts : argument describing the options. It can contains the fields:
%        -- args        : not considered if open_dialog is set to one. It
%                         must contains the fields:
%            --- iX       : indexes or names of the variable for which the
%                           sensitivity is computed. Any name not
%                           corresponding to a variable name or index out
%                           of the interval [1 P.DimX] will be ignored.
%            --- iP       : indexes or names of the parameter for which the
%                           sensitivity is computed. Names or indexes not
%                           corresponding to a parameter or an initial
%                           condition will be ignored.
%            --- tspan    : time instant(s) when the sensitivity is
%                           computed
%        -- props       : define the properties for which the sensitivity
%                         must be computed (DO NOT USE IT, THE
%                         FUNCTIONNALITY IS NOT USABLE)
%        -- taus        : column vector of size numel(props) x 1 or single
%                         value indicating the time instant at which the
%                         sensitivity of props must be computed.
%        -- open_dialog : (Optional, default=1) if set to one, open a
%                         dialog box asking for iP, iX and tspan. If used,
%                         the fields in args are only considered as
%                         suggestions for the dialog box.
%        -- plot_histo  : (Optional, default=1) if set to 1, plot the
%                         sensitivity
%        -- cutoff      : (Optional, default=0) cut off limit in percentage
%                         of the highest value
%        -- stat_type   : (Default='aver_sum') String defining the method
%                         used to compute the sensitivity. Can be either
%                         'aver_sum', in which case the provided
%                         sensitivity is the average of the sensitivity
%                         over all time point, and over the trajectories
%                         associated to the param set selected in ipts -
%                         the sensitivity of props in not computed (TO
%                         IMPLEMENT) ; or 'aver_end' in which case the
%                         answered sensitivity is the average of the
%                         sensitivity at time point tspan(end) over the
%                         trajectories associated to the parameter set
%                         selected in ipts. The sensitivity of the
%                         properties provided in opts.props is computed (TO
%                         UPDATE: THE FUNCTION QMITL_SEVALDIFF IS NOT UP TO
%                         DATE) ; or 'aver_max' ... TODO
%
% Outputs:
%   - M    : the values of sensitivities. The dimension of M is
%            numel(opts.iX)+numel(opts.props) x numel(opts.iP).
%   - opts : the structure of options used
%
% Example:
%
%
%See also
%

%
%  Plots histogram(s) of logarithmic sensitivities of state variables iX
%  w.r.t. parameters iP at a given time tspan for a parameter vector ipts
%  in P. iX, iP and tspan are provided through an input dialog box, except
%  when opts.args is given
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check inputs


if isempty(P.pts)
    error('SplotSensBar:NoPts','The field pts of P is empty !');
end

if ~exist('ipts','var')
    ipts=[];
end
if isempty(ipts)
    ipts = 1:size(P.pts,2);
end

% deal with options

% default values
open_dialog = 1;
cutoff = 0;
plots = 1;

if exist('opts','var')
    % do we need a dialog box to enter arguments (default: yes)
    if isfield(opts,'open_dialog')
        open_dialog = opts.open_dialog;
    end
        
    % cut off limit in percentage of the highest value (default: 0)
    if isfield(opts,'cutoff')
        cutoff = opts.cutoff;
    end
    
    if isfield(opts,'plot_histo')
        plots = opts.plot_histo;
    end
    
else
    opts = [];
end

if(open_dialog)
    % deal with dialog box for histogram axes
    try
        [iX, iP, tspan] = GetArgSensiBar(Sys.DimX, Sys.ParamList, opts.args);
    catch  %#ok<CTCH>
        try
            [iX, iP, tspan] = GetArgSensiBar(Sys.DimX, Sys.ParamList);
        catch ME
            warndlg(['Problem: ' ME.message ', giving up'] );
            close;
            return ;
        end
    end
    
    if isempty(tspan)
        if isfield(P, 'traj')
            tspan = P.traj(1).time;
        else
            error('SplotSensiBar:notspan','No precomputed trajectories.')
        end
    end
else
    iX = opts.args.iX;
    iP = opts.args.iP;
    tspan = opts.args.tspan;
end

% what type of computation do we do (default: average)
if isfield(opts,'stat_type')
    stat_type = opts.stat_type;
else
    if (numel(tspan)==1)
        stat_type = 'aver_end';  % compute sensitivities at a fixed time
        tspan = {0 tspan};
    else
        stat_type = 'aver_sum';
    end
end

% for properties evaluation
if isfield(opts, 'props')
    props = opts.props;
else
    props={};
end

if isfield(opts, 'taus')
    taus = opts.taus;
else
    taus = zeros(numel(props),1);
end

if(numel(taus)==1)
    taus = repmat(taus,numel(props),1);
end

if ~isnumeric(iX)
    iX = FindParam(P,iX);
    iX = iX(iX<=P.DimX); %keep only existing variables
end

if ~isnumeric(iP)
    iP = FindParam(P,iP);
    iP = iP(iP<=size(P.pts,1)); %keep only existing parameters
end

% From now on I shoud have Sys, ipts, tspan, iX, iP, prop, and taus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Recompute trajectories if needed

Ptmp = CreateParamSet(Sys, iP); 
Ptmp.pts = P.pts(:,ipts);
P = ComputeTrajSensi(Sys, Ptmp, tspan);

%  Compute the values for the histo bars

Mend = zeros(numel(iX)+numel(props), numel(iP));

switch (stat_type)
    case {'aver_sum'}  % integrate local sensitivities over the trajectories and average
        
        % 1/ Compute bars for variable sensitivities
        for i = 1:numel(ipts)
            
            traj = P.traj(i); %NM : shouldn't we use the traj_ref field? I don't think so, at least,
            time = traj.time;             % until we update ComputeTrajSensi
            
            for j = 1:numel(iX)
                for k = 1:numel(iP)
                    is = (find(P.dim==iP(k))-1)*size(traj.X,1)+iX(j);
                    
                    dx = traj.XS(is,:);  % dX/dp[t]
                    x = traj.X(iX(j),:);  % X[t]
                    
                    % replace zeros by small quantities
                    ind = find(abs(x)<1e-16);
                    x(ind) = sign(x(ind))*1e-16;
                    x(x==0) = 1e-16;
                    
                    p = traj.param(iP(k));    % p
                    xs = (dx*p)./abs(x);
                    
                    XS =  trapz(time, xs)/time(end); % computes the sum
                    
                    % Compute the sum
                    Mend(j,k) = Mend(j,k)+XS;
                end
            end
        end  % end i = ipts
        
        Mend = Mend/numel(ipts);
        
        % 2/ Compute bars for properties sensitivities
        
        % TODO
    
    %NM : maybe rename 'aver_end' in 'aver_time' or any name expressing
    % that the sensitivity is computed at a specific time
    case {'aver_end'}
        
        % 1/ Compute bars for variable sensitivities
        for i = 1:numel(ipts)
            traj = P.traj(i); %NM : shouldn't we use the traj_ref field? I don't think so, at least,
                                            % until we update ComputeTrajSensi
            for j = 1:numel(iX)
                for k = 1:numel(iP)
                    is = (find(P.dim==iP(k))-1)*size(traj.X,1)+iX(j);
                    
                    dx = traj.XS(is,end);  % dX/dp[t]
                    x = traj.X(iX(j),end);  % X[t]
                    
                    % replace zeros by small quantities
                    ind = find(abs(x)<1e-16);
                    x(ind) = sign(x(ind))*1e-16;
                    x(x==0) = 1e-16;
                    
                    p = traj.param(iP(k));    % p
                    xs = (dx*p)./abs(x);
                    
                    % Compute the average
                    Mend(j,k) = Mend(j,k)+xs;
                end
            end
        end  % end i = ipts
        
        Mend = Mend/numel(ipts);
        
        % 2/ Compute bars for properties sensitivities
        for j = 1:numel(props)
            for k = 1:numel(iP)
                %NM : should we compute the sensitivity at time tau or at
                %     the end of the simulation? I personally believe that
                %     we should compute the sensitivity at the end of the
                %     simulation.
                [p, x, dx] = QMITL_SEvalDiff(Sys, props{j}, P,  tspan, iP(k), taus(j)); %NM: use traj.time(end) instead of taus(j)?
                
                % replace zeros by small quantities
                ind = find(abs(x)<1e-16);
                x(ind) = sign(x(ind))*1e-16;
                x(x==0) = 1e-16;
                
                xs = (dx.*p)./abs(x);
                
                % Compute the average
                Mend(numel(iX)+j,k) = mean(xs);
                
            end
            
        end % end case aver_end
        
    case {'aver_max'}
        
        % 1/ Compute bars for variable sensitivities
        for i = 1:numel(ipts)
            traj = P.traj(i);
            
            for j = 1:numel(iX)
                for k = 1:numel(iP)
                    is = (find(P.dim==iP(k))-1)*size(traj.X,1)+iX(j);
                    %[dx idx] = max(abs(traj.XS(is,:)));
                    
                    dx = traj.XS(is, :);  % dX/dp[t]
                    x = traj.X(iX(j),:);  % X[t]
                    
                    % replace zeros by small quantities
                    ind = find(abs(x)<1e-16);
                    x(ind) = sign(x(ind))*1e-16;
                    x(x==0) = 1e-16;
                    
                    p = traj.param(iP(k));    % p
                    xs = (dx*p)./abs(x);
                    [~,idx] = max(abs(xs));
                    xs = xs(idx);
                    
                    % Compute the average
                    Mend(j,k) = Mend(j,k)+xs;
                end
            end
        end  % end i = ipts
        
        Mend = Mend/numel(ipts);
        
        % 2/ Compute bars for properties sensitivities
        for j = 1:numel(props)
            for k = 1:numel(iP)
                [p, x, dx] = QMITL_SEvalDiff(Sys, props{j}, P,  tspan, iP(k), taus(j)); %NM : shouldn't we compute all over tspan?
                
                % replace zeros by small quantities
                ind = find(abs(x)<1e-16);
                x(ind) = sign(x(ind))*1e-16;
                x(x==0) = 1e-16;
                
                xs = (dx.*p)./abs(x);
                
                % Compute the average
                Mend(numel(iX)+j,k) = mean(xs); %NM: the mean, really? not the max?
            end
        end
end % end switch

% Cut off negligible values
M = max(max(abs(Mend)));
Mend(abs(Mend)<cutoff*M) = 0;

if plots
    plot_histo(Mend,P, iX,props,iP);
end

end

function h = plot_histo3d(Mend, S, iX, props, iP)

h = figure;
h = bar3(Mend,0.5,'detached');

% Ticks labels

ytick_labels = cell(numel(iX)+numel(props));
for j = 1:numel(iX)
    ytick_labels(j) = S.ParamList(iX(j));
end
for j = 1:numel(props)
    ytick_labels(numel(iX)+j) = {get_id(props{j})};
end

xtick_labels = cell(1,numel(iP));
for k = 1:numel(iP)
    xlabel = S.ParamList{iP(k)};
    if (iP(k) <= S.DimX)
        xlabel = [xlabel '(0)'];
    end
    xtick_labels(k) = {xlabel};
end

set(gca, 'XTick', 1:size(Mend,2));
set(gca, 'YTick', 1:size(Mend,1));
set(gca, 'XTickLabel',  xtick_labels );
set(gca, 'YTickLabel',  ytick_labels );
axis([0 size(Mend,2)+1 0 size(Mend,1)+1]);

hx = get(gca, 'xlabel');
set(hx, 'Interpreter','none');
hy = get(gca, 'ylabel');
set(hy, 'Interpreter','none');

shading interp;
colormap cool;
%  colorbar;
for i = 1:length(h)
    zdata = get(h(i),'Zdata');
    set(h(i),'Cdata',zdata);
    set(h,'EdgeColor','k');
end

end

function h = plot_histo1(M, S, iX, props, iP)

h = figure;
%nb_histo = numel(iX)+numel(props);

% y labels

ytick_labels = cell(1,numel(iP));
for k = 1:numel(iP)
    ylabel = S.ParamList{iP(k)};
    if (iP(k)<= S.DimX)
        ylabel = [ylabel '(0)'];
    end
    ytick_labels(k) = {ylabel};
end

% plotting sensitivities of variables


barh(M);
set(gca, 'YTick', 1:numel(iP), 'YTickLabel',  ytick_labels);
hy = get(gca, 'ylabel');
set(hy, 'Interpreter','none');

legend(S.ParamList{iX});

end
