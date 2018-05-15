function [signals, monitors] = get_monitors(formulas)
% get_monitors determines stl monitors from a set of formulas, looking for templates alw, ev, until, alw A=>B,  etc

signals = {};
monitors = {};
for itfo = 1:numel(formulas)
    formula = formulas{itfo};
    
    if isa(formula, 'char')||isa(formula, 'STL_Formula')
        monitor = stl_monitor(formula);
    elseif isa(formula, 'stl_monitor')
        monitor = formula;
    end
    find_template();
    signals = [signals setdiff(monitor.signals_in, signals, 'stable')];
    monitors = [monitors {monitor}];
end


    function find_template()
        
        switch (get_type(monitor.formula))
            case  {'alw', 'always'}
                monitor = alw_monitor(formula);
            case {'ev', 'eventually'}
                monitor  = ev_monitor(formula);
        end
        
    end
end

