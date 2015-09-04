function PplotFormula(Sys, P, phis, ipts, break_level)
%PPLOTFORMULA plots formula evaluation. Each formula is plotted in a
% different subplot. The evaluation for all parameter vector is plotted in
% the same plot.
%
% Synopsis : PplotFormula(Sys, P, phis[, ipts[, break_level]])
% 
% Inputs:
%  - Sys         : the system
%  - P           : the parameter set. The evaluation of formulas described
%                  in phis must be computed.
%  - phis        : formulas which are plotted.
%  - ipts        : (Optional, default or empty=all parameter vectors)
%                  indexes of parameter vectors for which the evaluations
%                  of phis are plotted.
%  - break_level : (Optional,default=0) indicates the deepth until which
%                  formula described in phis are decomposed
% 
% Outputs:
%  - None, but a figure, hopefully, a nice one.
% 
%See also SEvalProp DiscrimPropValues
%


% manage inputs
if ~exist('break_level','var')
    break_level = 0;
end
if(break_level>0)
    phis_tmp = [];
    for i_pt = 1:numel(phis)
        broken_props = STL_Break(phis(i_pt),break_level);
        phis_tmp = [phis_tmp broken_props(:)]; %#ok<AGROW>
    end
    phis = phis_tmp;
end

if(~exist('ipts','var')||isempty(ipts))
    ipts = 1:size(P.pts,2);
end

% prepare figure
figure;
nb_phis = numel(phis);
if isfield(Sys,'time_mult')
    time_mult = Sys.time_mult;
else
    time_mult = 1;
end

% plot formula evaluation
for np = 1:numel(phis) % for each property
    phi = phis(np);
    i_phi = find_prop(P, get_id(phi));

    subplot(nb_phis, 1, np);
    hold on;
    xlabel('tau');
    title(disp(phi), 'Interpreter','none');
    
    for i_pt = ipts
        phi_tspan = P.props_values(i_phi,i_pt).tau;
        phi_val = P.props_values(i_phi,i_pt).val;
        plot(phi_tspan*time_mult, phi_val);
        plot([phi_tspan(1) phi_tspan(end)]*time_mult, [0 0],'-k');
        stairs(phi_tspan*time_mult, (phi_val>0)*max(abs(phi_val))/2,'-r','LineWidth', 4);
        YLim = get(gca, 'YLim');
        YLim(1) = min([-max(abs(phi_val))/2, YLim(1)]);
        set(gca,'YLim', YLim);
    end
    grid on;
    
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
