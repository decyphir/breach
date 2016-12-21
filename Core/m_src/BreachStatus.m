classdef BreachStatus < handle
    % BreachStatus a simple class to describe status of an object or output
    %              of an algorithm
    
    properties
        status = 0 % unique id for the status 
        statusList = {} % cell array of strings describing the status.
    end
    
    methods
        
        function addStatus(this, new_status, new_label)
            % Add new active status
            this.status = this.status + new_status;
            this.statusList = [this.statusList {new_label}];
        end
        
        function resetStatus(this)
            % Reset active status
            this.statusList = {};
        end
        
        function setStatus(this, new_status, new_label)
            % Set active status (reset previously active status)
            this.status = new_status;
            this.statusList = {new_label};
        end
        
        
        function st = getStatus(this)
            if isempty(this.statusList)
                st = 'No active status. \n';
            else
                n_status = numel(this.statusList);
                st = ['Active status = ' num2str(this.status) '\n'];
                for ist = 1:n_status
                    st = [st '- ' this.statusList{ist} '\n'];
                end
                
            end
        end
        
        
        function printStatus(this)
           fprintf(getStatus(this));
        end
        
    end
    
end