classdef BreachParam < handle
% BreachParam generic class to set and get parameters for simulation -
% mostly used to handle various workspaces, also structure and arrays 
    properties
        name
        ws='base'
    end
       
    methods 
        function this =BreachParam(name, ws)
            
            if nargin ==0 
                this.name= '';
                return;
            end
            
            this.name = name;
            if nargin>=2
                this.ws = ws;                
            end
        end
        
        function v = getValue(this)
            WS = this.getWorkspace();
            v = evalin(WS, this.name);
        end
        
        function setValue(this,v)
            WS = this.getWorkspace();
            assignin(WS, this.name,v);
        end
        
        function WS = getWorkspace(this)
            if exist(this.ws, 'file')==4   % this is a simulink model workspace
                warning('off', 'Simulink:Commands:SaveMdlWithDirtyWorkspace');
                WS = get_param(this.ws, 'ModelWorkspace');
            elseif isequal(this.ws, 'base')
                WS = 'base';
            else
                error(['Workspace ' this.ws ' invalid.']);
            end
        end
    end
end