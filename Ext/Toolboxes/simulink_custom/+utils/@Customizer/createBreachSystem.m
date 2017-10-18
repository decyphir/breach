% Copyright (c) 2008-2016 MonkeyProof Solutions B.V.
% Use is subject to the MIT license

function schema = createBreachSystem(callbackInfo) %#ok<INUSD> callbackInfo might be used in a later stadium
% Schema function: defines matchSizeOfBlocks menu item
 
schema              = sl_action_schema;     % Initialize schema 
schema.label        = 'Save and create a Breach interface';         % Set menu item label
schema.accelerator  = 'CTRL+ALT+B';         % Set accelerator/short-cut
schema.callback     = @createBreachSystemCb; % Set callback function 
 
end

function createBreachSystemCb(callbackInfo)
% Callback function: createBreachSystem menu item
 
mdl = gcs;
mdl = strsplit(mdl, '/');
mdl = mdl{1};
save_system(mdl);
Name = get_param(mdl, 'Name');
close_system(mdl);
B = evalin('base', ['BreachSimulinkSystem(''' Name ''')']);
signals = B.GetSignalList();
niou_sigs = select_cell_gui(signals, 'Select signals from list');
if isequal(niou_sigs,0) % pushed cancel
    open_system(mdl);
    return;
end

assignin('base','sigs__', niou_sigs);

params = B.GetParamsSysList();
if ~isempty(params)
    params = select_cell_gui(params, 'Select parameters from list');
    if isequal(params,0) % pushed cancel
        open_system(mdl);
        return;
    end
end
    assignin('base','params__', params);
    evalin('base', ['B=BreachSimulinkSystem(''' Name ''', params__,[],sigs__);']);
    evalin('base', 'B.RunGUI');
end