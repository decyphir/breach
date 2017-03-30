% Copyright (c) 2008-2016 MonkeyProof Solutions B.V.
% Use is subject to the MIT license

classdef Customizer < handle
    
    methods (Static)
        
        function childrenFcns = getCustomizeMethods(customizers)
            % Utility function to create childrenFcns out of the Customizers 
            % for custom menu schema
            
            childrenFcns    = {};
            
            % Loop over Customizers
            for iCustomizer = 1:length(customizers)
                
                % Retrieve class name for the customizer
                className           = class(customizers{iCustomizer});
                classNameEnd        = regexp(className, '\.(\w+)$', 'tokens', 'once');
                
                % Retrieve all methods implemented in the customizer, except for
                % constructor 'Customizer'
                methodNames         = setdiff(methods(customizers{iCustomizer}), [methods(mfilename); classNameEnd]);
                
                % Create function handles for each method
                customizerFunctions = cellfun(@(c)str2func([className '.' c]), methodNames, 'UniformOutput', false);
                
                % Add a separator between the function handles of each Customizer
                childrenFcns    	= [childrenFcns customizerFunctions.' 'separator']; %#ok<AGROW>
            end
            
            % Remove the last 'separator'
            childrenFcns(end)   = [];
            
        end
        
    end
    
end
