classdef BreachDomain
    % BreachDomain Implements types and fundamental sets
    %  Instead of having switch everywhere, should probably consider derived
    %  classes
  
    properties
        type='double'
        domain
        enum
    end
    
    methods
        function this = BreachDomain(type, domain)
        % BreachDomain(type, domain) 
        
            if nargin>1 && isempty(domain)
                domain=[];
            end
            
            switch nargin
                case 0
                    this.type  = 'double';
                    this.domain = [];
                case 1
                    if isnumeric(type)
                        if numel(type)==2
                            this.type = 'double';
                            this.domain = type;
                        else
                            error('BreachDomain:wrong_type', 'BreachDomain first argument should be the string ''double'', ''int'', ''enum'' or ''bool'', or an interval.')
                        end
                    else
                        this.type  = type;
                        this.domain = [];
                    end
                case 2
                    
                    if isnumeric(type)
                        if size(type,2)==2
                            this.type = 'double';
                            this.domain = type;
                        else
                            error('BreachDomain:wrong_type', 'BreachDomain first argument should be the string ''double'', ''int'', ''enum'' or ''bool'', or an interval.')
                        end
                    else
                        
                        this.type  = type;
                        this.domain = domain;
                    end
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
                case 'double'
                otherwise
                    error('BreachDomain:wrong_type', 'BreachDomain first argument should the string ''double'', ''int'', ''enum'' or ''bool'', or an interval.');
            end
        end
        
        function bool = is_default(this)
        % BreachDomain.is_default()
            bool = isequal(this.type, 'double') && isempty(this.domain);
        end
        
        function new_x  = checkin(this,x)
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
            if isempty(this.domain)
                new_x =  round(x);
            else
                new_x = min(max(round(x), ceil(this.domain(1))), floor(this.domain(2)));
            end
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
        
        function st = short_disp(this, not_if_default)
            % short disp
            if this.is_default()&& exist('not_if_default', 'var')&&not_if_default
                st= '';
                return;
            end
            
            if isempty(this.domain)
                st = ['of type ' this.type];
            else
                st = ['of type ' this.type ' in [' num2str(this.domain(1)) ', ' num2str(this.domain(2)) ']' ];
            end
            
        end
        
        function new_x = sample_grid(this, num_samples)
            new_x = linspace(this.domain(1),this.domain(2), num_samples);
            new_x = this.checkin(new_x);
        end
        
        
        function x = sample(varargin)
        % sample multi domain sampling. See BreachSet.SampleDomain
            
            %% process parameters
            
            % first arguments are domains
           i =2; 
            while isa(varargin{i}, 'BreachDomain')
                i = i+1;
            end
            domains = varargin(1:i-1);
            num_dom = i-1;
            
            % next ones are num_samples and then method
            num_samples = varargin{i};
            
            try
                method = varargin{i+1};
            end
            
            if ~exist('method','var')||isempty(method)
                method = 'rand';
            end
            
            % if all is selected, combine new samples
            combine_x=0;
            if iscell(num_samples)               
                for is =1:numel(num_samples)
                    combine_x = combine_x || isequal(num_samples{is}, 'all');  
                end
            end
            combine_x= combine_x||isequal(num_samples, 'all')||isequal(method,'grid')||isequal(method,'corners');
                
            if ischar(method)||isscalar(method)
                m = method;
                method = cell(1, num_dom);
                for ic = 1:num_dom
                    method{ic} = m;
                end
            end
            
            if ischar(num_samples)||isscalar(num_samples)
                ns = num_samples;
                num_samples = cell(1, num_dom);
                for ic = 1:num_dom
                    num_samples{ic} = ns;
                end
            end
            
            if isnumeric(num_samples)
                num_samples = num2cell(num_samples);
            end
            
            % creates new samples
            for ip = 1:numel(domains)
                dom = domains{ip};
                if isequal(num_samples{ip}, 'all')&&~isequal(method{ip}, 'corners')
                    x{ip} = dom.sample_all();
                else
                    if isequal(method{ip}, 'grid')
                        x{ip} = dom.sample_grid(num_samples{ip});
                    elseif isequal(method{ip}, 'rand')
                        x{ip} = dom.sample_rand(num_samples{ip});
                    elseif isequal(method{ip}, 'corners')
                        x{ip} = [dom.domain(1) dom.domain(2)]; 
                    end
                end
                num_x(ip) = numel(x{ip});
            end
            
            
            % Combine new samples
            if num_dom>1
                if combine_x
                    idx = N2Nn(num_dom, num_x);
                    for ip = 1:num_dom
                        new_x(ip,:) = x{ip}(1, idx(ip,:));
                    end
                    x = new_x;
                else
                    x = cell2mat(x');
                end
            else
                x= x{1};
            end
            
            
         end
        
        
    end
end