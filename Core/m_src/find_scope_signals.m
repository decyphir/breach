function signals = find_scope_signals(mdl)
% Limitations: 
%  - floating scopes signals are not captured by the current
%    mechanism. A workaround is to log signals vizualized in the floating
%    scope (in parameters of the scope, tab history)
%  - input signals for scopes need a name 

%
 
load_system(mdl);
scopes = find_system(mdl,'SearchDepth', inf, 'BlockType', 'Scope');

nb_scopes = numel(scopes);

signals={};
new_signals = {};

scope_save_names = {};

for i_scope = 1:nb_scopes
    scope = scopes{i_scope};
    scope_save_name = get_param(scope, 'SaveName');
    scope_save_name = genvarname(scope_save_name, scope_save_names);
    scope_save_names = {scope_save_names{:}, scope_save_name};
    set_param(scope, 'SaveName', scope_save_name);
    find_signals();
    signals = {signals{:} new_signals{:}};
end

% code for individual scopes, populate new_signals
    function find_signals()
        new_signals = {};
        in_signals = get_param(scope, 'InputSignalNames');
        lh = get_param(scope, 'LineHandles');
        lh = lh.Inport;
        
        for ins = 1:numel(in_signals)
            in_sig = in_signals{ins};
            
            % if input signals of the scope are unnamed, name them after
            % the scope name
            scope_name = get_param(scope,'Name');
            if isempty(in_sig)
                new_signal = genvarname(scope_name,{scope_name new_signals{:} signals{:}});
            else
                new_signal = in_sig;
            end
            % remove fancy characters from name
            new_signal = regexprep(new_signal,'\W','_');
            
            % ensures signal in has same name
            set(lh(ins), 'Name', new_signal);

            new_signals = {new_signal, new_signals{:}};
        end
        
        set_param(scope,'SaveToWorkspace', 'on');
        set_param(scope,'LimitDataPoints','off');
    end
end
