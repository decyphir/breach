function [mu, mustar, sigma, Pr] = SPropSensi(Sys, P, phi, opt)
% SPROPSENSI estimates the global sensitivity of a property wrt some
% parameters (Morris method). This function is adapted from the strategy
% described in "Global Sensitivity Analysis, A Primer", Saltelli et al,
% p113++
%
% Synopsis: [mu, mustar, sigma] = SPropSensi(Sys, P, phi, opt)
%
% Input:
%    - P  is a parameter set for Sys. P may contains many parameter sets,
%         but only the first will be considered, so it is recommanded that
%         P contains only one parameter set. The value of all parameters of
%         P, not in opt.params is defined by the first parameters values in
%         P.
%    - phi is an STL property
%    - opt is an option structure with the following fields :
%
%        - tspan       time domain computation of the trajectories
%        - tprop       time instant (scalar) when to eval prop satisfaction (default tspan(1))
%        - params      variable parameters
%        - lbound      lower bounds for the search domain
%        - ubound      upper bounds for the search domain
%        - p           number of levels (p-grid). Recommended to be even (default 4)
%        - r           number of trajectories (default 10)
%        - k           DEPRECATED ; REPLACED BY r
%        - plot        if 1, plots histograms (default 0)
%        - muGraphOpt      if plot=1, define graphical options for mu graph
%                          (optional)
%        - mustarGraphOpt  if plot=1, graphical options for mu* graph
%                          (optional)
%        - sigmaGraphOpt   if plot=1, graphical options for sigma graph
%                          (optional)
%
% Output:
%   - mu      expectation of elementary effects
%   - mustar  expectation of absolute values of elementary effects
%   - sigma   variance of elementary effects
%
% Example1 (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys);
%   phi = QMITL_Formula('phi','ev (x1[t] > 3)');
%   opt.tspan = 2:0.1:5;
%   opt.params = {'a','F'};
%   opt.lbound = [0.15, 5];
%   opt.ubound = [0.35, 25];
%   opt.p = 8;
%   opt.r = 100;
%   opt.plot = 1;
%   opt.muGraphOpt = {'XScale','log'};
%   [mu, mustar, sigma] = SPropSensi(Sys, P, phi, opt);
% 
% Example2 (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys);
%   P = SetParam(P, {'x1h', 'x1l', 'T'}, [.3, -.3, 5]);
%   [~,props] = QMITL_ReadFile('oscil_prop.stl');
%   phi = props{end};
%   opt.tspan = 0:0.1:5;
%   opt.params = {'a','F'};
%   opt.lbound = [0.15, 5];
%   opt.ubound = [0.35, 25];
%   opt.plot = 1;
%   [mu, mustar, sigma] = SPropSensi(Sys, P, phi, opt);
%
%See also QMITL_Formula QMITL_ReadFile CreateParamSet
%


if isfield(opt, 'tspan')
    tspan = opt.tspan;
elseif isfield(P, 'traj')
    tspan = P.traj(1).time;
elseif isfield(Sys, 'tspan')
    tspan = Sys.tspan;
else
    tspan = 0:.2:10;
end

if isfield(opt, 'tprop')
    tprop = opt.tprop;
else
    tprop = tspan(1);
end

if ~isfield(opt, 'p')
    opt.p = 4;
end

if ~isfield(opt, 'r')
    if isfield(opt,'k')
        opt.r = opt.k;
    else
        opt.r = 10;
    end
end

Sys.p = P.pts(:,1);
Sys.ParamList = P.ParamList;
Pr = CreateParamSet(Sys, opt.params, [opt.lbound' opt.ubound']);

Pr = pRefine(Pr, opt.p, opt.r);

%NM : it is better to compute the truth value of phi at time=tprop
if tprop < tspan(1)
    tspan = [tprop,tspan];
%elseif tprop > tspan(end)    % don't compute it in case it is too far from
%    tspan = [tspan, tprop];  % the last time instant computed
end

Pr = ComputeTraj(Sys, Pr, tspan);

[Pr, Y] = SEvalProp(Sys, Pr, phi, tprop);

[mu, mustar, sigma] = EEffects(Y, Pr.D, opt.p);


if opt.plot
    [~,isort] = sort(abs(mu));
    h = figure;
    subplot(3,1,1);
    barh(mu(isort));
    title('Expectation of elementary effects (mu)')
    set(gca, 'YTick', 1:numel(opt.params), 'YTickLabel', opt.params(isort));
    if isfield(opt,'muGraphOpt')
        set(gca,opt.muGraphOpt{:});
    end
    
    subplot(3,1,2);
    barh(mustar(isort));
    title('Expectation of absolute values of elementary effects (mu*)')
    set(gca, 'YTick', 1:numel(opt.params), 'YTickLabel', opt.params(isort));
    if isfield(opt,'mustarGraphOpt')
        set(gca,opt.mustarGraphOpt{:});
    end
    
    subplot(3,1,3);
    barh(sigma(isort));
    title('Standard deviation of elementary effects (sigma)')
    set(gca, 'YTick', 1:numel(opt.params), 'YTickLabel', opt.params(isort));
    if isfield(opt,'sigmaGraphOpt')
        set(gca,opt.sigmaGraphOpt{:});
    end
    
    fig_resize(h,1,2.5)
end

end
