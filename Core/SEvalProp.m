function [P, val] = SEvalProp(Sys, P, props, tau, ipts, bool_plot, break_level)
%SEVALPROP Eval property for previously computed trajectories
%
% Usage: [Pf, val] = SEvalProp(Sys, P, prop[ , tau[, ipts[, bool_plot[, break_level ]]]])
%
% Inputs:
%  - Sys         : system
%  - P           : Param set with trajectories
%  - prop        : QMITL property(ies)
%  - tau         : Time instant(s) when to estimate properties. If not
%                  provided, the time instants considered for computing the
%                  trajectory are used. It may be a scalar, in which case,
%                  all the formula are evaluated at this time point, or it
%                  may be an array of size 1 x numel(props), thus
%                  indicating the time point of evaluation of each formula.
%  - ipts        : Indices of parameter vectors for which the properties is
%                  evaluated (Optional, Default=all parameter sets).
%  - break_level : (Optional) defines the deep of breaking of props. If
%                  lower or equal to 1, it is ignored. If greater or equal
%                  to two, SEvalProp answers the evaluation of the props
%                  and all sub-formula of props until the deep provided.
%
% Outputs:
%  - Pf  : param set with prop_namse, prop and prop_values fields
%  - val : an array containing the quantitative satisfaction of properties
%          for each trajectory at the first time point of tau. The
%          dimension of val is numel(props) x numel(ipts)
%
%See also QMITL_Formula CreateParamSet ComputeTraj Sselect
%


% check arguments

if ~exist('break_level','var')
    break_level = 0;
end

if(break_level>0)
    nprops = [];
    for ii = 1:numel(props)
        nprops = [ nprops QMITL_Break(props(ii),break_level) ]; %#ok<AGROW>
    end
    props = nprops;
end

if ~exist('bool_plot','var')
    bool_plot = 0;
end

if(~exist('ipts','var')||isempty(ipts))
    ipts = 1:size(P.pts,2);
end

if ~exist('tau','var')
    tau = [];
elseif isscalar(tau)
    tau = ones(1,numel(props))*tau;
end

if ~iscell(props)
    props = {props};
end

if ~isfield(P,'props')
    P.props = [];
end

if ~isfield(P,'props_names')
    P.props_names = {} ;
end

if ~isfield(P,'traj_ref')
    P.traj_ref = 1:numel(P.traj);
end

% do things

% setup plots if needed

if(bool_plot)
    figure;
    nb_prop = numel(props);
    if isfield(Sys,'time_mult')
        time_mult = Sys.time_mult;
    else
        time_mult = 1;
    end
end

val = zeros(numel(props),numel(ipts)); %initialize array containing truth values for each property and each param set
props_values(1:numel(ipts)) = deal(struct()); % Temporary line containing the growing evaluation of the formula
for np = 1:numel(props) % for each property
    
    prop = props{np};  % prop = current property
    prop_name =  get_id(prop);
    iprop = find_prop(P,prop_name);
    
    if(bool_plot)
        subplot(nb_prop, 1, np);
        hold on;
        xlabel('tau');
        title(disp(prop), 'Interpreter','none');
    end
    
    if(~iprop)
        % if the property does not exist in P, we add it to P
        P.props_names = [P.props_names {prop_name}];
        P.props = [P.props prop];
        iprop = numel(P.props_names);
    end
    
    prop = QMITL_OptimizePredicates(Sys,prop);
    fprintf(['Checking ' prop_name  '\n'...
             '[             25%%           50%%            75%%               ]\n ']);
    iprog = 0; %idx of progression bar
    
    for ii = ipts % we compute the truch value of prop for each param set
        traj_tmp = P.traj(P.traj_ref(ii));
        Ptmp = Sselect(P,ii);
        if isempty(tau)
            tau_tmp = traj_tmp.time;
        else
            tau_tmp = tau(np);
        end
        props_values(ii).tau = tau_tmp;
        props_values(ii).val = QMITL_Eval(Sys, prop, Ptmp, traj_tmp, tau_tmp);
        val(np,ii) = props_values(iprop,ii).val(1);
        if isnan(val(np,ii))
            warning('SEvalProp:NaNEval','Warning: property evaluated to NaN');
        end
        
        while(floor(60*ii/numel(ipts))>iprog)
            fprintf('^');
            iprog = iprog+1;
        end
        
        % plot property values
        if(bool_plot)
            phi_tspan = props_values(ii).tau;
            phi_val = props_values(ii).val;
            plot(phi_tspan*time_mult, phi_val);
            plot([phi_tspan(1) phi_tspan(end)]*time_mult, [0 0],'-k');
            stairs(phi_tspan*time_mult, (phi_val>0)*max(abs(phi_val))/2,'-r');
            grid on;
        end
    end
    P.props_values(iprop,:) = props_values; % we copy all evaluation in once to avoid inconsistent parameter set
    
    fprintf('\n');
end

end

function idx = find_prop(P,st)
%FIND_PROP finds the index of a property in a parameter set.
%
% Synopsis: idx = find_prop(P,st)
%
% Input:
%  - P  the parameter set containing the evaluation of properties
%  - st a string describing the name of the searched property
%
% Output:
%  - the index of the property evaluation if found, 0 otherwize
%
try
    for idx = 1:numel(P.props_names)
        if strcmp(st,P.props_names{idx})
            return;
        end
    end
catch %#ok<CTCH>
end

idx=0;

end

