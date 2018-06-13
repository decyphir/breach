function [signals, monitors] = get_monitors(formulas)
% get_monitors determines stl monitors from a set of formulas, looking for templates alw, ev, until, alw A=>B,  etc

signals = {};
monitors = {};
for itfo = 1:numel(formulas)
    formula = formulas{itfo};
    if ischar(formula)
        formula = STL_Formula(STL_NewID('req'), formula);
    end
    
    if isa(formula, 'STL_Formula')
        % checks whether we have parameter constraint
        sigs  = STL_ExtractSignals(formula);       
        if isempty(sigs)
            monitor = param_constraint_monitor(formula);
        else
            monitor = stl_monitor(formula);
        end
    elseif isa(formula, 'stl_monitor')
        monitor = formula;
    end
    if ~isa(monitor, 'param_constraint_monitor')
        find_template();
        signals = union(signals, monitor.signals, 'stable');
    end
    monitors = [monitors {monitor}];
end

    function find_template()
        p0 = monitor.p0;
        switch (get_type(monitor.formula))
            case  {'alw', 'always'}
                child = get_children(monitor.formula);
                if strcmp(get_type(child{1}), '=>')
                    monitor = alw_A_implies_B_monitor(formula);
                end
                monitor = alw_monitor(formula);
            case {'ev', 'eventually'}
                monitor  = ev_monitor(formula);
            case {'until'}
            otherwise  % default to top alw if horizon is 0
               hor = get_horizon(formula);
               if hor==0
                   alw_formula = STL_Formula(['alw_' get_id(monitor.formula)], ['alw ' get_id(monitor.formula)]);
                   monitor = alw_monitor(alw_formula);
                   monitor.p0 = p0;
               end
        end
    end
end

