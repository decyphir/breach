classdef BreachDomain
    % Instead of having switch everywhere, should probably consider derived
    % classes
    % checkin does not take into account enum and tabu_list yet ..
    properties
        type='double'
        domain
        enum
    end
    
    methods
        function this = BreachDomain(type, domain)
            if nargin>1 && isempty(domain)
                domain=[];
            end
            
            switch nargin
                case 0
                    this.type  = 'double';
                    this.domain = [];
                case 1
                    this.type  = type;
                    this.domain = [];
                case 2
                    this.type  = type;
                    this.domain = domain;
             end
            
            % init enum (I feel this will explode on me some day)
            switch this.type
                case 'int'
                    if numel(this.domain)~=2
                       this.enum = round(this.domain);
                       this.domain  = [min(this.enum), max(this.enum)];
                    else
                        this.enum = this.domain(1):this.domain(2);
                    end
                case 'enum'
                    this.enum = this.domain;
                    this.domain = [min(this.enum), max(this.enum)];
                case 'bool'
                    this.domain = [0 1];
                    this.enum = [0 1];
            end
        end
        
        function new_x = checkin(this,x)
            switch this.type
                case 'int'
                    new_x = this.checkin_int(x);
                case 'bool'
                    new_x = this.checkin_bool(x);
                case 'enum'
                    new_x = this.checkin_enum(x);
                    
                case 'double'
                    new_x =   this.checkin_double(x); % mostly out of bounds
            end
        end
        
        function new_x = checkin_int(this,x)
            new_x = min(max(round(x), ceil(this.domain(1))), floor(this.domain(2)));
        end
        
        function x = checkin_double(this,x)
            if ~(isempty(this.domain))
                x = min(max(x,this.domain(1)), this.domain(2));
            end
        end
        
        function x = checkin_bool(this,x)
            x = (x~=0);
        end
        
        function x = checkin_enum(this,x)            
            for ix= 1:numel(x)
                [~ , imin] = min(abs(this.enum-x(ix))); 
                x(ix) = this.enum(imin);
            end
          end
        
        
        function all_x = sample_all(this)
            all_x = this.enum;
        end
        
        function new_x = sample_rand(this, num_samples)
        % assumes bounded domain 
            switch this.type
                case {'enum', 'int'}
                    new_x = this.enum(randi(numel(this.enum),1,num_samples));
                case 'double'
                    new_x = (this.domain(2)-this.domain(1))*rand(1,num_samples) + this.domain(1);
            end
        end
        
        function new_x = sample_grid(this, num_samples)
            new_x = linspace(this.domain(1),this.domain(2), num_samples);
            new_x = this.checkin(new_x);
        
        end
        
    end
end