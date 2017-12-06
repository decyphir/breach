classdef BreachArrayParam < BreachParam
    properties
        idx =[1 1]
    end
    
    methods
        
        
        function this = BreachArrayParam(name, idx, ws)
            if nargin==0
                return;
            end
            this.name = name;
            this.idx = idx;
            switch nargin
                case 3
                    this.ws = ws;
            end
        end
        
        function setValue(this,v)
            WS = this.getWorkspace();
            ar =evalin(WS,  this.name);
            ar(this.idx(1),this.idx(2)) = v;
            assignin(WS, this.name, ar);
        end
        
        function  v= getValue(this)
            WS = this.getWorkspace();
            ar =evalin(WS,  this.name);
            v = ar(this.idx(1),this.idx(2));
        end
    end
end
