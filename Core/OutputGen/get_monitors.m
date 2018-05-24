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
        monitor = stl_monitor(formula);
    elseif isa(formula, 'stl_monitor')
        monitor = formula;
    end
    find_template();
    signals = [signals setdiff(monitor.signals_in, signals, 'stable')];
    monitors = [monitors {monitor}];
end

    function find_template()
        p0 = monitor.p0;
        
        switch (get_type(monitor.formula))
            case  {'alw', 'always'}
                monitor = alw_monitor(formula);
            case {'ev', 'eventually'}
                monitor  = ev_monitor(formula);
            case {'until'}
            otherwise  % default to top alw if horizon is 0
               hor = get_horizon(formula);
               if hor==0
                   alw_formula = STL_Formula(['alw_' get_id(monitor.formula)], ['alw ' get_id(monitor.formula)]);
                   monitor = alw_monitor(alw_formula);
               end
        end
    end
end

