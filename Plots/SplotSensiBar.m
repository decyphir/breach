function [Mend, opts] =  SplotSensiBar(Sys, P, ipts, opts)
%SPLOTSENSIBAR plots an histogram of local sensitivities
% 
% Synopsis:  [M, opts] = SplotSensiBar(Sys, P, ipts [, opts])
% 
% Inputs:
%  - Sys  : the system
%  - P    : a parameter set. It must contain a field pts or an error is
%           thrown. It may contain many parameter vectors.
%  - ipts : (Optional, default or empty=all parameter vectors) array
%           providing the indexes of the parameter vectors of P to consider
%           for computing the sensitivity. The sensitivity is averaged over
%           all these parameter vectors.
%  - opts : argument describing the options. It can contains the fields:
%     -- plot_histo  : (Optional, default=1) if set to 1, open a figure to
%                      plot the sensitivity
%     -- open_dialog : (Optional, default=1) if set to one, open a dialog
%                      box asking for iP, iX and tspan. If used, the fields
%                      in args are only considered as suggestions for the
%                      dialog box.
%     -- args        : only considered as suggestion if open_dialog is set
%                      to one. It can contain the fields:
%         --- iX       : (Optional, default=all) indexes or names of the
%                        variable for which the sensitivity is computed.
%                        Any name which does not match a variable name or
%                        any index not in the interval [1, P.DimX] will be
%                        ignored.
%         --- iP       : (Optional, default=P.dim) indexes or names of the
%                        parameter for which the sensitivity is computed.
%                        Names not matching a parameter or variable name
%                        will be ignored, as well as indexes not in [1,
%                        size(P.pts,1)].
%         --- tspan    : time point(s) for trajectories computation. May be
%                        not provided if Sys contains a field tspan.
%                        Otherwise, an error will *very likely* be thrown.
%     -- phis        : (default=empty) array of STL_Formula defining the
%                      properties for which the sensitivity must be
%                      computed.
%     -- props       : DEPRECATED, use phis instead.
%     -- taus        : column vector of size numel(phis) x 1 or single
%                      value indicating the time instant at which the
%                      sensitivity of props must be computed.
%     -- cutoff      : (Optional, default=0) cut off limit in percentage of
%                      the highest value: all sensitivity lower than cutoff
%                      times the highest sensitivity will be set to 0.
%     -- stat_type   : (Default='aver_sum') String defining the method
%                      used to compute the sensitivity. Can be either
%                      'aver_sum', 'aver_time' or 'aver_max'. If it is
%                      'aver_sum', the sensitivity to parameters is the
%                      average over the trajectory and then, the average
%                      over all trajectories. The sensitivity to formula is
%                      set to zero because it is not implemented (TO
%                      IMPLEMENT). If stat_type is 'aver_end' or
%                      'aver_time', the function computes the sensitivity
%                      of variables at the time point tspan(end), averaged
%                      over all the trajectories. The sensitivity of
%                      formulas is computed at time point provided by taus
%                      and averaged over all trajectories. If stat_type is
%                      'aver_max', for each trajectory, the sensitivity
%                      with the highest absolute value over the trajectory
%                      is keept and all sensitivity are average over all
%                      trajectories. For STL formula, the sensitivity is
%                      computed at each time instant in tspan and the
%                      highest is keept. Then, the sensitivity is averaged
%                      over all trajectories.
%
% Outputs:
%  - M    : the values of sensitivities. The dimension of M is
%           numel(opts.iX)+numel(opts.props) x numel(opts.iP).
%  - opts : the structure of options used
% 
% Example (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys, {'x1','G'}, [-10 10; 0 5],4);
%   P = SetParam(P, {'x1h', 'x1l', 'T'}, [.3, -.3, 5]);
%   [~,phis] = STL_ReadFile('oscil_prop.stl');
%   phis = [phis{[1,end]}] % keep the first and the last formula
%   opts.plot_histo = 1;
%   opts.open_dialog = 0;
%   opts.args.iX = {'x2'};
%   opts.args.iP = {'F','x1'};
%   opts.args.tspan = 0:0.1:10;
%   opts.phis = phis;
%   opts.taus = 0;
%   opts.stat_type = 'aver_time';
%   Sensi = SplotSensiBar(Sys,P,[],opts)
% 
%See also PphiSensiLocal
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
    error('SplotSensBar:emptyPtsField','The field pts of P is empty !');
end

if(~exist('ipts','var') || isempty(ipts))
    ipts = 1:size(P.pts,2);
end

% deal with options

if ~exist('opts','var')
    opts = [];
end

if isfield(opts,'plot_histo')
    plots = opts.plot_histo;
else
    plots = 1; %default is to plot sensi
end

% do we need a dialog box to enter arguments (default: yes)
if isfield(opts,'open_dialog')
    open_dialog = opts.open_dialog;
else
    open_dialog = 1;
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
            return;
        end
    end
    
    if isempty(tspan)
        if isfield(P, 'traj')
            tspan = P.traj{1}.time; %NM: not correct because two trajectories may have different time field
        else
            error('SplotSensiBar:noTspan','No precomputed trajectories.')
        end
    end
else
    if ~isfield(opts.args,'iX')
        iX = 1:P.DimP;
    else
        iX = opts.args.iX;
    end
    if isfield(opts.args,'iP')
        iP = opts.args.iP;
    else
        iP = P.dim;
    end
    if isfield(opts.args,'tspan')
        tspan = opts.args.tspan;
    elseif isfield(Sys,'tspan')
        tspan = Sys.tspan;
    end
end
if ~isnumeric(iX)
    iX = FindParam(P,iX);
end
iX = iX(iX<=P.DimX); %keep only existing variables
iX = iX(iX>0);
if ~isnumeric(iP)
    iP = FindParam(P,iP);
end
iP = iP(iP<=size(P.pts,1)); %keep only existing parameters
iP = iP(iP>0);

% for properties evaluation
if isfield(opts, 'phis')
    phis = opts.phis;
elseif isfield(opts, 'props')
    warning('SplotSensiBar:deprecatedUseOfProps',...
        'The use of the field opts.props is deprecated. Please, use opts.phis.');
    phis = opts.props;
else
    phis = [];
end
if iscell(phis)
    phis = [phis{:}];
end

if isfield(opts, 'taus')
    taus = opts.taus;
    if(numel(taus)==1)
        taus = repmat(taus,numel(phis),1);
    end
else
    taus = zeros(numel(phis),1);
end

% cut off limit in percentage of the highest value (default: 0)
if isfield(opts,'cutoff')
    cutoff = opts.cutoff;
else
    cutoff = 0;
end

% what type of computation do we do (default: average)
if isfield(opts,'stat_type')
    stat_type = opts.stat_type;
else
    if(numel(tspan)==1)
        stat_type = 'aver_end';  % compute sensitivities at a fixed time
        tspan = {0 tspan};
    else
        stat_type = 'aver_sum';
    end
end

% From now on we shoud have Sys, ipts, tspan, iX, iP, phis, and taus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Recompute trajectories if needed

P = Sselect(P, ipts);
old_dim = P.dim;
P = SAddUncertainParam(P, iP);
P = SDelUncertainParam(P, old_dim, iP);
if ~isempty(iX)
    P = ComputeTrajSensi(Sys, P, tspan);
else
    P = ComputeTraj(Sys, P, tspan);
end

%  Compute the values for the histo bars

Mend = zeros(numel(iX)+numel(phis), numel(iP));

switch(stat_type)
    case 'aver_sum'  % integrate local sensitivities over the trajectories and average
        
        % 1/ Compute variable sensitivities
        for ii = 1:size(P.pts,2)
            traj = P.traj{P.traj_ref(ii});
            time = traj.time;
            
            for jj = 1:numel(iX)
                for kk = 1:numel(iP)
                    is = (find(P.dim==iP(kk))-1)*size(traj.X,1)+iX(jj);
                    
                    dx = traj.XS(is,:);  % dX/dp[t]
                    x = traj.X(iX(jj),:);  % X[t]
                    
                    % replace zeros by small quantities
                    ind = find(abs(x)<1e-16);
                    x(ind) = sign(x(ind))*1e-16;
                    x(x==0) = 1e-16;
                    
                    p = traj.param(iP(kk));    % p
                    xs = (dx*p)./abs(x);
                    
                    XS = trapz(time, xs)/time(end); % computes the average over the trajectory
                    
                    % Compute the sum
                    Mend(jj,kk) = Mend(jj,kk)+XS;
                end
            end
        end  % end ii = P.pts
        
        Mend = Mend/size(P.pts,2); % average over all trajectories
        
        % 2/ Compute formula sensitivities
        Mend(numel(iX)+1:end,:) = PPhiSensiLocal(Sys,P,phis,tspan,taus,iP,ipts,stat_type,cutoff);
    
    case {'aver_time','aver_end'}
        
        % 1/ Compute variable sensitivities
        for ii = 1:size(P.pts,2)
            traj = P.traj{P.traj_ref(ii});
            for jj = 1:numel(iX)
                for kk = 1:numel(iP)
                    is = (find(P.dim==iP(kk))-1)*size(traj.X,1)+iX(jj);
                    
                    dx = traj.XS(is,end);  % dX/dp[t]
                    x = traj.X(iX(jj),end);  % X[t]
                    
                    % replace zeros by small quantities
                    ind = find(abs(x)<1e-16);
                    x(ind) = sign(x(ind))*1e-16;
                    x(x==0) = 1e-16;
                    
                    p = traj.param(iP(kk));    % p
                    xs = (dx*p)./abs(x);
                    
                    % Sum all to compute the average
                    Mend(jj,kk) = Mend(jj,kk)+xs; % sum of sensitivities at the end of the trajectories
                end
            end
        end  % end ii = P.pts
        
        Mend = Mend/size(P.pts,2); % average over all trajectories
        
        % 2/ Compute formula sensitivities
        Mend(numel(iX)+1:end,:) = PPhiSensiLocal(Sys,P,phis,tspan,taus,iP,ipts,'aver_time',cutoff);
%         for jj = 1:numel(phis)
%             
%             [p, x, dx] = STL_SEvalDiff(Sys, phis(jj), P, tspan, iP, taus(jj));
%             
%             % replace zeros by small quantities
%             ind = abs(x)<1e-16;
%             x(ind) = sign(x(ind))*1e-16;
%             x(x==0) = 1e-16;
%             x = repmat(x,numel(iP),1);
%             
%             xs = (dx.*p)./abs(x);
%             
%             Mend(numel(iX)+jj,:) = mean(xs,2)'; % Compute the average over all trajectories
%             
% %             for kk = 1:numel(iP) % for each parameter
% %                 [p, x, dx] = STL_SEvalDiff(Sys, phis(jj), P, tspan, iP(kk), taus(jj));
% %                 
% %                 % replace zeros by small quantities
% %                 ind = abs(x)<1e-16;
% %                 x(ind) = sign(x(ind))*1e-16;
% %                 x(x==0) = 1e-16;
% %                 
% %                 xs = (dx.*p)./abs(x);
% %                 
% %                 Mend(numel(iX)+jj,kk) = mean(xs); % Compute the average over all trajectories
% %             end
%             
%         end % end jj = phis
        
    case 'aver_max'
        
        % 1/ Compute variable sensitivities
        for ii = 1:size(P.pts,2)
            traj = P.traj{P.traj_ref(ii});
            for jj = 1:numel(iX)
                for kk = 1:numel(iP)
                    is = (find(P.dim==iP(kk))-1)*size(traj.X,1)+iX(jj);
                    %[dx idx] = max(abs(traj.XS(is,:)));
                    
                    dx = traj.XS(is, :);  % dX/dp[t]
                    x = traj.X(iX(jj),:);  % X[t]
                    
                    % replace zeros by small quantities
                    ind = find(abs(x)<1e-16);
                    x(ind) = sign(x(ind))*1e-16;
                    x(x==0) = 1e-16;
                    
                    p = traj.param(iP(kk));    % p
                    xs = (dx*p)./abs(x);
                    [~,idx] = max(abs(xs));
                    xs = xs(idx); % keep the maximal sensitivity
                    
                    Mend(jj,kk) = Mend(jj,kk)+xs; % Compute the average
                end
            end
        end  % end ii = P.pts
        
        Mend = Mend/size(P.pts,2);
        
        % 2/ Compute formula sensitivities
        Mend(numel(iX)+1:end,:) = PPhiSensiLocal(Sys,P,phis,tspan,taus,iP,ipts,stat_type,cutoff);
%         for jj = 1:numel(phis)
%             for kk = 1:numel(iP)
%                 xs_max = zeros(1,size(P.pts,2)); %zero is the lowest absolute value
%                 for tau = tspan
%                     [p, x, dx] = STL_SEvalDiff(Sys, phis(jj), P, tspan, iP(kk), tau);
%                     
%                     % replace zeros by small quantities
%                     ind = find(abs(x)<1e-16);
%                     x(ind) = sign(x(ind))*1e-16;
%                     x(x==0) = 1e-16;
%                     
%                     xs = (dx.*p)./abs(x);
%                     
%                     idx = abs(xs_max)<abs(xs);
%                     xs_max(idx) = xs(idx); % keep the max
%                 end
%                 Mend(numel(iX)+jj,kk) = mean(xs); % average over all trajectories
%             end
%         end % end jj = phis
end % end switch

% Cut off negligible values
M = max(max(abs(Mend)));
Mend(abs(Mend)<cutoff*M) = 0;

if plots
    plot_histo(Mend, P, iX, phis, iP, taus);
end

end
