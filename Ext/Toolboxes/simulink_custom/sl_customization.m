function sl_customization(cm)
% Register custom menu in the Simulink Editor's menu bar.
cm.addCustomMenuFcn('Simulink:MenuBar', @getCustomSimulinkMenu);
end


function schemaFcns = getCustomSimulinkMenu(callbackInfo)
% Custom menu function: returns schema functions
schemaFcns = {@customMenu};
end


function schema = customMenu(callbackInfo) 
% Schema function: defines the custom menu

% Initialize schema
schema = sl_container_schema; 

% Set menu label
schema.label  = 'Breach'; 

% Initialize Customizers to add
customizers = {utils.Customizer()}; % , styleguide.Customizer()};

% Generate childrenFcns for schema
schema.childrenFcns = Customizer.getCustomizeMethods(customizers);

end