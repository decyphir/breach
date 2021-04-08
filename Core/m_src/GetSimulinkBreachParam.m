function bparam = GetSimulinkBreachParam(mdl, pname)
% GetBreachParam returns a BreachParam object accessor for setting and getting value of parameter pname for the simulator

spl = strsplit(pname, '.');
if numel(spl) == 1
    % scalar or array
    [bb, array_pname, idx] = isarray();
    if ~bb
        [ws] = getWorkspace(pname);
        bparam = BreachParam(pname, ws);
    else
        [ws] = getWorkspace(array_pname);
        bparam= BreachArrayParam(array_pname, idx, ws);
    end
else % we have a struct
    struct_pname = spl{1};
    field_pname = pname(numel(struct_pname)+2:end);
    [ws] = getWorkspace(struct_pname);
    [~, ~, idx] = isarray();
    bparam = BreachStructParam(struct_pname, field_pname,idx, ws);
    
end

    function [bb, array_pname, idx] = isarray()
        tok = regexp(pname,'^(\w+)__at_(\d+)_(\d+)$', 'tokens');
        array_pname = '';
        bb = ~isempty(tok);
        idx  = [1 1];
        if bb
            array_pname = tok{1}{1};
            idx(1) = str2double(tok{1}{2});
            idx(2) = str2double(tok{1}{3});
        end
    end

    function [ws] = getWorkspace(vname)
        % look into base workspace
        var = evalin('base', ['exist(''' vname ''',''var'')']);
        if var~=0
            ws = 'base';
        else            
            var= Simulink.findVars(mdl, 'WorkspaceType','model', 'Name', vname);
            if ~isempty(var)
                ws = mdl;
            else
                error('BreachSimulinkSystem:param_not_found', ['Parameter ' vname ' is not defined in the model nor base workspace.' ]);                
            end                        
        end
    end
end
