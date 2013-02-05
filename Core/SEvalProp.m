function [S, val] = SEvalProp(Sys, S, props, tau, ipts, bool_plot, break_level)
%
%   SEVALPROP Eval property for previously computed trajectories
%
%   Usage: [Pf, val] = SEvalProp(Sys, Ptraj, prop [ , tau, ipts ])
%
%   Inputs:
%
%    - Sys         system
%    - Ptraj       param set with trajectories
%    - prop        property(ies)
%    - tau         time instant(s) when to estimate properties. If not
%                  provided, the time instants considered for computing the
%                  trajectory are used.
%    - ipts        indices of trajectories for which to eval properties
%
%
%   Outputs:
%
%    - Pf          param set with prop_values field
%    - val         an array containing the quantitative satisfaction of
%                    properties for each trajectory
%


% check arguments

if (~exist('ipts','var')||isempty(ipts))
    ipts = 1:size(S.pts,2);
end

if (~exist('break_level','var'))
    break_level = 0;
end

if (break_level>0)
    nprops = [];
    for i = 1:numel(props)
        nprops = [ nprops QMITL_Break(props(i),break_level) ];
    end
    props = nprops;
end

if ~isfield(S,'props')
    S.props = [];
end

if ~isfield(S,'props_names')
    S.props_names = {} ;
end

if ~isfield(S,'traj_ref')
    S.traj_ref = 1:numel(S.traj);
end

if (~exist('tau','var')||isempty(tau))
    tau0 = [];
else
    tau0 = tau;
end

if (~exist('bool_plot','var'))
    bool_plot = 0;
end

% do things

% setup plots if needed

if (bool_plot)
    figure;
    nb_prop = numel(props);
    if (isfield(Sys,'time_mult'))
        time_mult = Sys.time_mult;
    else
        time_mult = 1;
    end
end

val = zeros(numel(props),numel(ipts)); %initialize array containing truth values for each property and each param set
for np = 1:numel(props) % for each property
    
    prop = props(np);  % prop = current property
    prop_name =  get_id(prop);
    iprop = find_prop(S,prop_name);
    
    if (bool_plot)
        subplot(nb_prop, 1, np);
        hold on;
        xlabel('tau');
        title(disp(prop), 'Interpreter','none');
    end
    
    if ~iprop
        % if the property does not exist in S, we add it to S
        S.props_names = [S.props_names {prop_name}];
        S.props = [S.props prop];
        iprop = numel(S.props_names);
    end
    
    prop = QMITL_OptimizePredicates(Sys,prop);
    fprintf(['Checking ' prop_name  '\n'...
             '[             25%%           50%%            75%%               ]\n ']);
    iprog = 0; %idx of progression bar
    
    Ptmp = Sselect(S,1); % copie S en ne gardant que le premier parameter set
    
    for i = ipts % we compute the truch value of prop for each param set
        while (floor(60*i/numel(ipts))>iprog)
            fprintf('^');
            iprog = iprog+1;
        end
        
        traj = S.traj(S.traj_ref(i));
        Ptmp.pts = S.pts(:,i); % we copy in Ptmp the ith param set ; BETTER TO USE Ptmp = Sselect(S,i) ???
        if isempty(tau0)
            tau = traj.time; % no need of "else tau=tau0" because tau0 is a copy of tau
        end
        S.props_values(iprop,i).tau = tau;
        S.props_values(iprop,i).val = QMITL_Eval(Sys,prop,Ptmp, traj, tau);
        val(np,i) = S.props_values(iprop,i).val(1);
        
        % plot property values
        if (bool_plot)
            phi_tspan = S.props_values(iprop,i).tau;
            phi_val = S.props_values(iprop,i).val;
            plot(phi_tspan*time_mult, phi_val);
            plot([phi_tspan(1) phi_tspan(end)]*time_mult, [0 0],'-k');
            stairs(phi_tspan*time_mult, (phi_val>0)*max(abs(phi_val))/2,'-r');
            grid on;
        end
        
    end
    
    fprintf('\n');
end


function i = find_prop(S,st)

i = 0;
for k = 1:numel(S.props_names)
    if strcmp(st,S.props_names{k})
        i = k;
        return;
    end
end
