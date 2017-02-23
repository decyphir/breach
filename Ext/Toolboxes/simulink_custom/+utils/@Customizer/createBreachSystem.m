% Copyright (c) 2008-2016 MonkeyProof Solutions B.V.
% Use is subject to the MIT license

function schema = createBreachSystem(callbackInfo) %#ok<INUSD> callbackInfo might be used in a later stadium
% Schema function: defines matchSizeOfBlocks menu item
 
schema              = sl_action_schema;     % Initialize schema 
schema.label        = 'Creates a Breach interface';         % Set menu item label
schema.accelerator  = 'CTRL+ALT+B';         % Set accelerator/short-cut
schema.callback     = @createBreachSystemCb; % Set callback function 
 
end

function createBreachSystemCb(callbackInfo)
% Callback function: createBreachSystem menu item
 
mdl = gcs;
Name = get_param(mdl, 'Name');
B = evalin('base', ['BreachSimulinkSystem(''' Name ''')']);
signals = B.GetSignalList();
niou_sigs = select_cell_gui(signals, 'Select signals from list');
assignin('base','sigs__', niou_sigs);

params = B.GetParamsSysList();
niou_params = select_cell_gui(params, 'Select parameters from list');
assignin('base','params__', niou_params);
evalin('base', ['B=BreachSimulinkSystem(''' Name ''', params__,[],sigs__);']);
evalin('base', 'B.RunGUI');
end