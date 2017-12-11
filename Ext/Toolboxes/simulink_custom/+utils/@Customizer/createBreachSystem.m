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
 
mdl = bdroot;
save_system(mdl);

BreachSimulinkWizard(mdl, [], 'RunFromSimulink', true);
end