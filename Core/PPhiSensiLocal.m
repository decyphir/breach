function [M, val] = PPhiSensiLocal(Sys, P, phis, tspan, taus, params, ipts, stat_type, cutoff, VERBOSE)
%PPHISENSILOCAl computes local sensitivity of STL formulas
% 
% Synopsis:  [M, val] = PPhiSensiLocal(Sys, P, phis, tspan[, taus[, params[, ipts[, stat_type[, cutoff[, VERBOSE]]]]]])
% 
% Inputs:
%  - Sys       : the system
%  - P         : a parameter set. It must contain a field pts or an error
%                is thrown. It may contain many parameter vectors. The
%                trajectories does not need to be computed, nor the
%                evaluation of phis.
%  - phis      : array of QMITL_Formula defining the properties for which
%                the sensitivity is computed.
%  - tspan     : time point(s) for trajectories computation. May be not
%                provided (aka: set to empty) if Sys contains a field
%                tspan.
%  - taus      : (Optional, default or empty=first time instant of tspan)
%                column vector of size numel(phis) x 1 or single value
%                indicating the time instant at which the sensitivity of
%                phis must be computed. Not used if stat_type is set to
%                'aver_max'
%   - params   : (Optional, default or empty=P.dim) indexes or names of the
%                parameter for which the sensitivity is computed. Names not
%                matching a parameter or variable name will be ignored, as
%                well as indexes not in [1, size(P.pts,1)].
%  - ipts      : (Optional, default or empty=all parameter vectors) array
%                providing the indexes of the parameter vectors of P to
%                consider for computing the sensitivity. The sensitivity is
%                averaged over all these parameter vectors.
%  - stat_type : (Optional, default or empty='aver_time') String defining
%                the method used to compute the sensitivity. Can be either
%                'aver_sum', 'aver_time' or 'aver_max'. If stat_type is
%                'aver_time', the sensitivity of formulas is computed at
%                time point provided by taus and averaged over all
%                trajectories. If stat_type is 'aver_max', the formula
%                sensitivity is computed at each time instant in tspan and
%                the highest is keept. Then, the sensitivity is averaged
%                over all trajectories. The 'aver_sum' case is not
%                implemented.
%  - cutoff    : (Optional, default=0) cut off limit in percentage of
%                the highest value: all sensitivity lower than cutoff
%                times the highest sensitivity will be set to 0.
%  - VERBOSE   : (Optional, default=true) Boolean indicating if the
%                computation progress bar is shown.
% 
% Outputs:
%  - M   : the values of sensitivities. The dimension of M is
%          numel(phis) x numel(params).
%  - val : the considered evaluation of phis at time taus. The size of val
%          is numel(phis) x numel(ipts)
% 
% Example (Lorentz84):
%   CreateSystem
%   P = CreateParamSet(Sys, {'x1','G'}, [-10 10; 0 5],4);
%   [~,phis] = QMITL_ReadFile('oscil_prop.stl');
%   phis = [phis{[1,end]}] % keep the first and the last formula
%   Sensi = PPhiSensiLocal(Sys, P, phis, 0:0.1:10, 0, {'F','x1'})
%   
%   P2 = CreateParamSet(Sys, 'x1', [1, 1]);
%   P2 = SetParam(P2, 'x1h', .3);
%   phi1 = phis(1) % == x1_high == x1[t]>x1h
%   [M, rho] = PPhiSensiLocal(Sys, P2, phi1, 0:0.1:10, 0, 'x1')
%   P2 = SetParam(P2, 'x1', 1.01);  % x1 is 1% more
%   P2 = ComputeTraj(Sys,P2,0:0.1:10);
%   rho*(1+M/100)  % expected rho
%   [~,rho] = SEvalProp(Sys, P2, phi1, 0) % got  1.4286% more :)
% 
%See also SPropSensi SplotSensiBar
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check inputs
if isempty(P.pts)
    error('SplotSensBar:emptyPtsField','The field pts of P is empty !');
end

if iscell(phis) % manage cell array (in case of...)
    phis = [phis{:}];
end

if(~exist('tspan','var')||isempty(tspan))
    if isfield(Sys,'tspan')
        tspan = Sys.tspan;
    else
        error('PPhiSensiLocal:noTspan','tspan is defined neither as parameter nor is Sys.');
    end
end

if(~exist('taus','var')||isempty(taus))
    taus = tspan(1);
end
if(numel(taus)==1)
    taus = repmat(taus,numel(phis),1);
end

if(~exist('params','var')||isempty(params))
    params = P.dim;
end
if ~isnumeric(params)
    params = FindParam(P,params);
end
params = params(params<=size(P.pts,1)); %keep only existing parameters
params = params(params>0);

if(~exist('ipts','var') || isempty(ipts))
    ipts = 1:size(P.pts,2);
else
    ipts = ipts(ipts>0);
    ipts = ipts(ipts<=size(P.pts,2));
end

% what type of computation do we do (default: average)
if(~exist('stat_type','var') || (~strcmp(stat_type,'aver_max')&&~strcmp(stat_type,'aver_sum')) )
    stat_type = 'aver_time';
end

% cut off limit in percentage of the highest value (default: 0)
if ~exist('cutoff','var')
    cutoff = 0;
end

% do we show progress bar
if ~exist('VERBOSE','var')
    VERBOSE = true;
end

% From now on we shoud have Sys, ipts, tspan, iX, iP, phis, and taus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = Sselect(P, ipts);

M = zeros(numel(phis), numel(params));
val = zeros(numel(phis), numel(ipts));

switch(stat_type)
    case 'aver_sum'
        % TODO
        fprintf('\nSorry, this is still in developpment...\n');
    
    case 'aver_time'
        for ii = 1:numel(phis)
            [p, x, dx] = QMITL_SEvalDiff(Sys, phis(ii), P, tspan, params, taus(ii), VERBOSE);
            
            % replace zeros by small quantities
            ind = abs(x)<1e-16;
            x(ind) = sign(x(ind))*1e-16;
            x(x==0) = 1e-16;
            
            val(ii,:) = x; % store evaluation of phis(ii) for all parameter vectors
            
            x = repmat(x,numel(params),1);
            xs = (dx.*p)./abs(x);
            
            M(ii,:) = mean(xs,2)'; % Compute the average over all trajectories
        end
        
    case 'aver_max'
        for ii = 1:numel(phis)
            for jj = 1:numel(params)       %TODO: CAN BE OPTIMIZED
                xs_max = zeros(1,size(P.pts,2)); % zero is the lowest absolute value
                for tau = tspan
                    [p, x, dx] = QMITL_SEvalDiff(Sys, phis(ii), P, tspan, params(jj), tau, VERBOSE);

                    % replace zeros by small quantities
                    ind = find(abs(x)<1e-16);
                    x(ind) = sign(x(ind))*1e-16;
                    x(x==0) = 1e-16;

                    xs = (dx.*p)./abs(x);
                    
                    idx = abs(xs_max)<abs(xs);
                    xs_max(idx) = xs(idx); % keep the max
                    
                    val(ii,ind) = x(ind); % keep the formula evaluation corresponding to the maximal sensitivity
                end
                M(ii,jj) = mean(xs_max); % average over all trajectories
            end
        end
end

% Cut off negligible values
Mmax = max(max(abs(M)));
M(abs(M)<cutoff*Mmax) = 0;

end
