classdef BreachDomain
    % BreachDomain Implements types and fundamental sets. 
    %
    %  TODO error handling when intersection of domain and enum is empty 
    % 
     
    properties
        type='double'     % can be 'int', 'bool', 'enum', 'double'
        domain            % always be an interval, empty means singleton
        enum              % only  used with type enum and 'bool'
    end
    
    methods
        function this = BreachDomain(type, domain,enum)
            % BreachDomain(type, domain, enum)
            %
            
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
                    elseif ismember(type, {'double', 'int' ,'bool'} )
                        this.type  = type;
                        this.domain = [];
                    elseif isequal(type, 'enum')
                        error('BreachDomain:enum_needs_values', 'Domain of type enum requires that values are provided in a second or third argument.');
                    else
                        error('BreachDomain:wrong_args', 'BreachDomain wrong arguments: unknown type or invalid domain.');
                    end
                    
                case 2
                    switch type
                        case {'double', 'bool'}
                            this.type = type;
                            if isempty(domain)||(size(domain,2)==2)
                                this.domain = domain;
                            else
                                error('BreachDomain:wrong_type', 'BreachDomain second argument should be empty or an interval.')
                            end
                        case 'int'
                            this.type = 'int';
                            if  all([0 2] ~= numel(domain))
                                this = BreachDomain('enum', [], domain);
                            else
                                this.domain = domain;
                            end
                        case {'enum'}
                             this.type = type;
                             if isnumeric(domain)&& (size(domain, 1) ==1|| size(domain, 2) ==1)
                                if size(domain, 2) ==1
                                    domain = domain';
                                end
                                this.enum = domain;
                                this.domain = [min(domain) max(domain)];
                             else
                                 error('BreachDomain:wrong_domain_or_enum', 'BreachDomain wrong arguments: domain or enum size invalid.');
                             end
                        otherwise
                            error('BreachDomain:wrong_args', 'BreachDomain wrong arguments: unknown type or invalid domain.');
                    end
                case 3 % only used for enum: define a dom
                    
                    switch type 
                        case 'enum'  % TODO some additional consistency checking
                            this.type = 'enum';
                            this.domain = domain;
                            this.enum = enum; 
                        otherwise
                            error('BreachDomain:wrong_args', 'BreachDomain wrong arguments: unknown type or invalid domain or enum.');
                    end
            end
            
            % init enum for bool 
           if isequal(this.type, 'bool')
                    this.enum = [0 1];
           end
        end
        
        function bool = is_default(this)
        % BreachDomain.is_default()
            bool = isequal(this.type, 'double') && isempty(this.domain);
        end
        
        function [new_x, is_in, dist]   = checkin(this,x)
        % BreachDomain.checkin(x) checks      
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
            
            if nargout >= 2
                    is_in = isequal(x,new_x);
            end
            if nargout>= 3               
                dist = abs(x-new_x);                    
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
            assert(~isempty(this.domain), ...
                    'BreachDomain:sample_all:empty_domain',...
                    'sample_all cannot sample an empty domain.');
            
          switch this.type
                case 'bool'
                    all_x = [false true];
                case 'int'
                    all_x = this.domain(1):this.domain(2);
                case 'enum'
                    all_x = this.enum;
                case 'double'
                    error ('BreachDomain:sample_all:double_domain',...
                    'sample_all cannot sample a double domain.');
          end
        end
            
          
        
                     
        function new_x = sample_rand(this, num_samples)
            % assumes bounded domain
            assert(~isempty(this.domain), ...
                    'BreachDomain:sample_rand:empty_domain',...
                    'sample_rand cannot sample an empty domain.');
            
            switch this.type
                case 'bool'
                    new_x = boolean(randi([0 1], 1, num_samples));
                case 'int'
                    new_x = randi([this.domain(1) this.domain(2)],1, num_samples);
                case 'enum'
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
                if ismember(this.type,{'bool', 'enum'})
                    st = [ st ' const in {' num2str(this.enum(1:min(10,end)))];
                    if numel(this.enum)>10
                        st = [st ' ... ' num2str(this.enum(end)) '}'];
                    else
                        st = [st '}'];
                    end
                    
                end
                
            else
                st = ['of type ' this.type ' in [' num2str(this.domain(1)) ', ' num2str(this.domain(2)) ']' ];
                if ismember(this.type,{'bool', 'enum'})
                    st = sprintf([st ' intersect with {'  num2str(this.enum(1:min(10,end)))]);
                    if numel(this.enum)>10
                        st = [st ' ... ' num2str(this.enum(end)) '}'];
                    else
                        st = [st '}'];
                    end
                    
                end
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
            
            if numel(varargin)>i
                i = i+1;
                method = varargin{i};
            end
            
            if ~exist('method','var')||isempty(method)
                method = 'rand';
            end
            
            % Last one is max_num_samples
            if numel(varargin)>i
                i = i+1;
                max_num_samples = varargin{i};
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
                    else
                        error('Unknown sampling method: %s', method{ip});
                    end
                end
                num_x(ip) = numel(x{ip});
            end
            
            
            % Combine new samples  --> this is where we need some smarts 
            if num_dom>1
                if combine_x
                    if  exist('max_num_samples','var')
                        idx = N2Nn(num_dom, num_x, max_num_samples);
                    else
                        idx = N2Nn(num_dom, num_x);
                    end
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