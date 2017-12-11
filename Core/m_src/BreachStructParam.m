classdef BreachStructParam < BreachArrayParam
    properties
        struct_name
        field_name
    end
    
    methods
        function this = BreachStructParam(struct_name,field_name,idx, ws)
            this.struct_name = struct_name;
            this.field_name = field_name;
            this.name = [struct_name '.' field_name];
            switch nargin
                case 2
                    this.idx = [1 1];
                case 3
                    this.idx=idx;
                case 4
                    this.idx = idx;
                    this.ws = ws;
            end
        end
        
        function setValue(this,v)
            WS = this.getWorkspace();
            st =evalin(WS,  this.struct_name);
            st.(this.field_name)(this.idx(1), this.idx(2)) = v;
            assignin(WS, this.struct_name, st);
        end
        
        function  v= getValue(this)
            WS = this.getWorkspace();
            st =evalin(WS,  this.struct_name);
            v= st.(this.field_name)(this.idx(1), this.idx(2));
        end
        
    end
end
