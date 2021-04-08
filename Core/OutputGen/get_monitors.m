function [signals, monitors] = get_monitors(formulas)
% get_monitors determines stl monitors from a set of formulas, looking for templates alw, ev, until, alw A=>B,  etc

signals = {};
monitors = {};
if ~iscell(formulas)
    formulas = {formulas};
end

for itfo = 1:numel(formulas)
    formula = formulas{itfo};
    if ischar(formula)
        formula = STL_Formula(STL_NewID('req'), formula);
    end
    if isa(formula, 'req_monitor')
        monitor = formula;
    else
        monitor = stl_monitor(formula);
    end
    
    %   if isa(formula, 'STL_Formula')
    %       % checks whether we have parameter constraint   % ONHOLD
    %        sigs  = STL_ExtractSignals(formula);
    %         if isempty(sigs)  %
    %             monitor = param_constraint_monitor(formula);
    %         else
    %           monitor = stl_monitor(formula);
    %         end
    %     elseif isa(formula, 'stl_monitor')
    %        monitor = formula;
    %   end
    %   if ~isa(monitor, 'param_constraint_monitor')
    if isa(monitor,'stl_monitor')
        find_template();
    end
    signals = union(signals, monitor.signals, 'stable');
    %   end
    monitors = [monitors {monitor}];
end

    function find_template()
        p0 = monitor.p0;
        switch (get_type(monitor.formula))
            case  {'alw', 'always'}
                if ~isa(monitor, 'alw_monitor')
                    monitor = alw_monitor(formula);
                end
%                 child = get_children(monitor.formula);
%                 if strcmp(get_type(child{1}), '=>')
%                     monitor = alw_A_implies_B_monitor(formula);
%                 else
%                     monitor = alw_monitor(formula);
%                 end
%            case {'ev', 'eventually'}
%                monitor  = ev_monitor(formula);
%             case {'until'}
%             otherwise  % default to top alw if horizon is 0
%                hor = get_horizon(formula);
%                if hor==0
%                    alw_formula = STL_Formula(['alw_' get_id(monitor.formula)], ['alw ' get_id(monitor.formula)]);
%                    if strcmp(get_type(formula), '=>')
%                        monitor = alw_A_implies_B_monitor(alw_formula); 
%                    else
%                        monitor = alw_monitor(alw_formula);
%                    end
%                end
        end
        monitor.p0 = p0;
    end
end

