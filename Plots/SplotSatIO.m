function [P, val] =  SplotSatIO(Sys, P, phis, depth, tau, ipts, inout, relabs)
%SplotSatIO computes and plots the satisfaction function wrt time of
% properties for computed trajectories
% 
% Synopsis: [Pf, val] = SplotSatIO(Sys, Ptraj, phis[, depth[, tau[, ipts]]])
% 
% Inputs:
%  - Sys   : the system
%  - Ptraj : parameter set with trajectories
%  - phis  : property(ies)
%  - depth : (Optional, default=0) computes and plot satisfaction of
%            subformulas up to depth
%  - tau   : (Optional, default or empty=traj.time for each trajectories)
%            time instant(s) when to estimate all properties
%  - ipts  : (Optional, default=all parameter vectors) parameter vectors
%            for which trajectories are considered to eval properties
%  - inout : (Optional, default='out') subset of variables over which the
%            quantitative semantics is computed: 'in' or 'out'.
%  - relabs : (Optional, default='rel') default robustness value for 
%            variables over which only a qualitative value is computed:
%            'rel' uses +inf/-inf yielding the relative robustness;
%            'abs' uses +0/-0 yielding the absolute robustness.
% 
% Outputs:
%  - Pf  : param set with prop_values field
%  - val : quantitative satisfaction of properties
%

InitBreach;
global BreachGlobOpt


% check arguments
if(~exist('ipts','var')||isempty(ipts))
    ipts = 1:size(P.pts,2);
end

if(~exist('tau','var')||isempty(tau)) % also manage empty cell
    tau = [];
end

if ~exist('depth','var')
    depth = inf;
end

if ischar(phis)
    phi_tmp__ = STL_Formula('phi_tmp__', phis );
    [P,val] = SplotSat(Sys,P, phi_tmp__ , depth, tau, ipts);
    return;
end

if(depth>0)
    nphis = [];
    for ii = 1:numel(phis)
        nphis = [nphis STL_Break(phis(ii), depth) ]; %#ok<*AGROW>
    end
    phis = nphis;
end

if ~isfield(P,'props')
    P.props = [];
    npb = 0;
else
    npb = numel(P.props);
end
if ~isfield(P,'props_names')
    P.props_names = {} ;
end
if ~isfield(P,'traj_ref')
    P.traj_ref = 1:numel(P.traj);
end

if ~exist('inout','var')
    inout = 'out';
end

if ~exist('relabs','var')
    relabs = 'rel';
end

%% setup plots if needed
nb_phis = numel(phis);
if isfield(Sys,'time_mult')
    time_mult = Sys.time_mult;
else
    time_mult = 1;
end

Ptmp = Sselect(P,1); % this avoid to use Sselect npb*ipts times
for np = npb+1:nb_phis+npb
    
    phi = phis(np-npb);
    phi_name = get_id(phi);
    iphi = find_prop(P, phi_name);
    if(~iphi)
        P.props_names = [P.props_names, phi_name];
        P.props = [P.props, phi];
        iphi = numel(P.props_names);
    end
    
    subplot(nb_phis, 1, np-npb);
    hold on;
    title(disp(phi), 'Interpreter', 'none');
    
    phi = STL_OptimizePredicates(Sys, phi);
%   fprintf(['Checking ' phi_name  '\n[             25%%           50%%            75%%               ]\n ']);
%   iprog = 0;
    
    for ii = ipts
        traj = P.traj{P.traj_ref(ii)};
        Ptmp.pts = P.pts(:,ii);
        if ~isempty(tau)
            phi_tspan = tau;
            phi_val = STL_Eval_IO(Sys, phi, Ptmp, traj, inout, relabs, tau);
            val(np-npb,ii) = phi_val(1);
        else
            [phi_val, phi_tspan] = STL_Eval_IO(Sys, phi, Ptmp, traj, inout, relabs);
            val(np-npb,ii) = phi_val(1);
        end
        
%        while(floor(60*ii/numel(ipts))>iprog)
%           fprintf('^');
%            iprog = iprog+1;
%        end
        %display(phi);
        %display(phi_tspan);
        %display(phi_val);
        
        % output property values
        P.props_values(iphi,ii).tau = phi_tspan;
        P.props_values(iphi,ii).val = phi_val;
        
        %plot(phi_tspan*time_mult, phi_val);
        %stairs(phi_tspan*time_mult, (phi_val>0)*max(abs(phi_val))/2,'-r');
        tsc = phi_tspan*time_mult;
        plot_style = 'stairs';
        
        if isfield(BreachGlobOpt, 'disable_robust_linear_interpolation')
            if BreachGlobOpt.disable_robust_linear_interpolation==0
                plot_style = 'plot';
            end
        end
        
        phi_bool = NaN(size(phi_val));
        for i = 1:size(phi_val,1)
            for j = 1:size(phi_val,2)           
                % Notice that the case phi_val(i,j) == 0 should yield 
                % phi_bool(i,j) == NaN. This is important for consistency 
                % of results in the case of absolute robustness
                if phi_val(i,j)>0 % numerical instability ??
                    phi_bool(i,j) = 1;
                elseif phi_val(i,j)<0 % numerical instability ??
                    phi_bool(i,j) = 0;
                end
            end
        end

        ax = plotyy(tsc, phi_val, tsc, phi_bool, plot_style, 'stairs' );
        set(ax(2), 'YLim', [-0.1 1.1], 'YTick', [0 1], 'YTickLabel', {'false', 'true'});
        if np-npb == 1
            legend('Quant. sat', 'Bool. sat');
        end
        grid on;
        
    end
    
%    fprintf('\n');
end
xlabel('time');
ax = get(gcf, 'Children');
ax = ax(arrayfun(@(c)(isa(c,'matlab.graphics.axis.Axes')), ax)==1);
linkaxes(ax, 'x');
h = zoom;
set(h,'Motion','horizontal','Enable','off');
            

end

function idx = find_prop(P, st)
%FIND_PROP finds the index of a property in a parameter set.
%
% Synopsis: idx = find_prop(P,st)
%
% Input:
%  - P  the parameter set containing the evaluation of properties
%  - st a string describing the name of the searched property
%
% Output:
%  - the index of the property evaluation if found, 0 otherwise
%

if ~isfield(P,'props_names')
    idx = 0;
    return;
else
    for idx = 1:numel(P.props_names)
        if strcmp(st,P.props_names{idx})
            return;
        end
    end
end

idx=0; % in case it is not found

end
