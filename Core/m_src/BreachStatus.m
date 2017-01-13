classdef BreachStatus < handle
    % BreachStatus a simple class to describe status of an object or output
    %              of an algorithm
    
    properties
        status = [] % unique integer id for the status 
        statusMap = containers.Map() % map with key and msg strings describing the status.
    end
    
    methods
        
        function addStatus(this, status, key, msg)
            % Add new active status
            this.status = [this.status status];
            
            if (~ischar(key) && numel(key)>20)
               error('key should be a string of max. 20 characters.'); 
            end
            
            this.statusMap(key) = msg;
        end
        
        function resetStatus(this)
            % Reset active status
            this.status = [];
            this.statusMap= containers.Map();
        end
        
        function setStatus(this, status, key, msg)
            % Set active status (reset previously active status)
            this.status = status;
            this.statusMap(key) = msg;
        end
        
        function st = getStatus(this)
            if isempty(this.statusMap)
                st = 'No active status. \n';
            else
                n_status = numel(this.statusMap.keys);
                keys= this.statusMap.keys;
                st = ['Active status = ' num2str(this.status) '\n'];
                for ist = 1:n_status
                    new_st = ['[' keys{ist} ']'];
                    msg = this.statusMap(keys{ist});
                    new_st(24:24+numel(msg)-1) = msg;
                    st = [st new_st  '\n'];
                end
                
            end
        end
        
        function printStatus(this)
           fprintf(getStatus(this));
        end
        
    end
    
end