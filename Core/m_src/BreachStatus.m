classdef BreachStatus < handle
    % BreachStatus a simple class to describe status of an object or output
    %              of an algorithm
    
    properties
        status = [] % unique integer id for the status 
        statusMap = containers.Map() % map with key and msg strings describing the status.
        verbose = 1 
        logged_msg
        max_logged_msg=100
    end
    
    methods
        
        function  my_name = whoamI(this)
            S = evalin('base', 'who');
            my_name= '__Nobody__'; 
            for iv = 1:numel(S)
              if isequal(evalin('base', S{iv}), this)
                  my_name = S{iv};
              end
            end
        end
        
        function new = copy(this)
            % copy operator, works with R2010b or newer.
            objByteArray = getByteStreamFromArray(this);
            new = getArrayFromByteStream(objByteArray);
        end

        
        function addStatus(this, status, key, msg)
            % Add new active status
            this.status = [this.status status];
            
            if (~ischar(key) && numel(key)>200)
               error('key should be a string of max. 200 characters.'); 
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
      
        function disp_msg(this, msg, verbose_min) 
            % disp_msg display a message if this.verbose is more than
            % second arg
            if nargin<3
                verbose_min = 1;
            end
            
            if (this.verbose>=verbose_min)
               disp(msg);                 
            end
            this.logged_msg = [this.logged_msg {msg}];
            
        end
        
    end
    
end