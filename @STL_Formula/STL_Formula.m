function [phi, phistruct] = STL_Formula(varargin)
%QSTL_FORMULA STL formula constructor
%
%     STL_Formula(id) creates a new formula named id
%
%     | STL_Formula(id,'phi_expression') where phi obeys the following grammar
%
%      phi_expression := predicate_expression | unary_op phi_expression | '('phi_expression')' binary_op '('phi_expression')'
%
%      predicate_expression :=   matlab_expression = matlab_expression '|' threshold = num_value, max_true_value = num_value
%                              | matlab_expression comp matlab_expression
%
%      comp := > | >= | < | <=
%
%      unary_op := not | unary_op_temp | unary_op_temp_[ti,tf]
%
%      unary_op_temp :=  ev | alw | eventually | always
%
%      binary_op  := or | and | until | until_[matlab_expession, matlab_expression]
%
%
%     | STL_Formula(id,'unary_op',phi) where 'unary_op' is one of 'not',
%                'eventually' or 'always' (abreviations 'ev' or 'alw')
%
%     | STL_Formula(id,'unary_op',phi, interval) where 'unary_op' is
%                either 'eventually' or 'always'(abreviations 'ev' or 'alw')
%
%     | STL_Formula(id,'binary_op',phi1, phi2) where 'binary_op' is
%                either 'or','and' or 'until'
%
%     | STL_Formula(id,'until',phi1,interval, phi2)
%
evalin('base','InitBreach');

global BreachGlobOpt

% test if formula already exists
if(nargin==2)
    % here we copy a formula 
    if isa(varargin{2},'STL_Formula')
        phi.id = varargin{1};
        phistruct = struct(phi);
        BreachGlobOpt.STLDB(phi.id) = phi;
        return;
    elseif ischar(varargin{2}) % here we reference an existing formula 
    st = varargin{2};
    if isKey(BreachGlobOpt.STLDB,st)
        phi = BreachGlobOpt.STLDB(st);
        return;
    end
    else 
       error('STL_Formula:Bad_argument_type', 'Second argument should be a string or a formula.');
    end
end

% OK, new formula or erroneous args, let's proceed
if(numel(varargin)==0)
    varargin{1} = 'phi';
end

phi.id = varargin{1};
phi.st = '';
phi.evalfn = [];
phi.interval = [0 inf];
phi.phi = [];
phi.phi1 = [];
phi.phi2 = [];
phi.phin = [];
phi.type = '';
phi.params = struct;
phi.params_interval = struct;
varargin = varargin(2:end);

switch numel(varargin)
    case 0
        phi.st = 'true';
        phi.evalfn = @true_formula;
        phi.type = 'predicate';
        
    otherwise
        st = varargin{1};
        if(numel(varargin)==1 && ischar(st))
            st = regexprep(st,'eventually', 'ev');
            st = regexprep(st,'always','alw');
            st = regexprep(st,'<=', '<');
            st = regexprep(st,'>=', '>');
            % deals with true and false
            st = regexprep(st, 'true', 'inf>0');
            st = regexprep(st, 'false', 'inf<0');
            varargin{1}= st;
        end
        
        phi = STL_Parse(phi,varargin{:});
end


phi=class(phi, 'STL_Formula');
phi = check_params(phi);
phistruct = struct(phi);
BreachGlobOpt.STLDB(phi.id) = phi;

end

function phi = STL_Parse(phi,varargin)
%STL_PARSE fills the field type, phi, phi1 and phi2
%
% Synopsis: phi = STL_Parse(phi, phi_str)
%
%      OR : phi = STL_Parse(phi, unary_op, phi0)
%      OR : phi = STL_Parse(phi, unary_op2, interv, phi2)
%      OR : phi = STL_Parse(phi, binary_op, phi1, phi2)
%      OR : phi = STL_Parse(phi, 'until', phi1, interv, phi2)
%      OR : phi = STL_Parse(phi, 'andn', [phi1, phi2, ..., phin])
%
% Inputs:
%  - phi       : is the STL_Formula to create
%  - phi_str   : a string describing the formula. This string must follow
%                the grammar described in the STL_Formula documentation
%  - phi0      : a STL Formula
%  - phi1      : a STL Formula
%  - phi2      : a STL Formula
%  - phin      : a STL Formula
%  - unary_op  : is either 'not', 'ev', 'alw', 'eventually' or 'always'
%  - unary_op2 : is either 'ev', 'alw', 'eventually' or 'always'
%  - binary_op : is either 'or', 'and' or '=>'
%  - interv    : is an interval
%
% Output:
%  - phi : a STL Formula structure
%

global BreachGlobOpt

switch(numel(varargin))
    
    case 1 % Here, the formula is defined by a string. We parse this string
        if ~ischar(varargin{1})
            error('STL_Formula:STL_Parse','Invalid formula');
        end
        st = varargin{1};
        st = regexprep(st, '^\s*', '');
        
        %% test or
        [success, st1, st2] = parenthesisly_balanced_split(st, '\<or\>');
        if success
            phi1 = STL_Formula([phi.id '1__'],st1);
            phi2 = STL_Formula([phi.id '2__'],st2);
            STLDB_Remove([phi.id '1__']);
            STLDB_Remove([phi.id '2__']);
            phi = STL_Parse(phi,'or', phi1, phi2);
            return
        end
        
        %% test implies
        [success, st1, st2] = parenthesisly_balanced_split(st, '\<=>\>');
        if success
            phi1 = STL_Formula([phi.id '1__'],st1);
            phi2 = STL_Formula([phi.id '2__'],st2);
            STLDB_Remove([phi.id '1__']);
            STLDB_Remove([phi.id '2__']);
            phi = STL_Parse(phi,'=>', phi1, phi2);
            return
        end
        
        %% test and
        [success, st1, st2] = parenthesisly_balanced_split(st, '\<and\>');
        if success
            phi1 = STL_Formula([phi.id '1__'],st1);
            phi2 = STL_Formula([phi.id '2__'],st2);
            STLDB_Remove([phi.id '1__']);
            STLDB_Remove([phi.id '2__']);
            
            phi = STL_Parse(phi,'and', phi1, phi2);
            return
        end
        
        %% test until
        [success, st1, st2] = parenthesisly_balanced_split(st, '\<until\>');
        interval = '[0 inf]';
        if success
            phi1 = STL_Formula([phi.id '1__'],st1);
            phi2 = STL_Formula([phi.id '2__'],st2);
            STLDB_Remove([phi.id '1__']);
            STLDB_Remove([phi.id '2__']);
            
            phi = STL_Parse(phi,'until', phi1, interval, phi2);
            return
        end
        
        %% test until_[ti,tf]
        [success, st1, st2, interval] = parenthesisly_balanced_split(st, '\<until_\[(.+?)\]\>');
        if success
            phi1 = STL_Formula([phi.id '1__'],st1);
            phi2 = STL_Formula([phi.id '2__'],st2);
            STLDB_Remove([phi.id '1__']);
            STLDB_Remove([phi.id '2__']);
            
            phi = STL_Parse(phi,'until', phi1, interval, phi2);
            return
        end
        
        %% test eventually
        [success, st1, st2] = parenthesisly_balanced_split(st, '\<ev\>');
        if success && isempty(st1)
            phi1 = STL_Formula([phi.id '1__'],st2);
            STLDB_Remove([phi.id '1__']);
            phi = STL_Parse(phi, 'ev', phi1);
            return
        end
        
        %% test eventually_[ti,tf]
        [success, st1, st2, interval] = parenthesisly_balanced_split(st, '\<ev_\[(.+?)\]\>');
        if success && isempty(st1)
            phi1 = STL_Formula([phi.id '1__'],st2);
            STLDB_Remove([phi.id '1__']);
            phi = STL_Parse(phi,'ev',interval,phi1);
            return
        end
        
        %% test always
        [success,st1, st2] = parenthesisly_balanced_split(st, '\<alw\>');
        if success && isempty(st1)
            phi1 = STL_Formula([phi.id '1__'],st2);
            STLDB_Remove([phi.id '1__']);
            phi = STL_Parse(phi, 'alw', phi1);
            return
        end
        
        %% test alw_[ti,tf]
        [success, st1, st2, interval] = parenthesisly_balanced_split(st, '\<alw_\[(.+?)\]\>');
        if success && isempty(st1)
            phi1 = STL_Formula([phi.id '1__'],st2);
            STLDB_Remove([phi.id '1__']);
            phi = STL_Parse(phi,'alw',interval,phi1);
            return
        end
        
        
        %% test not
        [success,st1, st2] = parenthesisly_balanced_split(st, '\<not\>');
        if success && isempty(st1)
            phi1 = STL_Formula([phi.id '1__'],st2);
            STLDB_Remove([phi.id '1__']);
            phi = STL_Parse(phi, 'not', phi1);
            return
        end
        
        % test expr op expr | params
        
        % parse additional params
        tokens = regexp(st, '(.+)\s*\|\s*(.+)','tokens');
        if ~isempty(tokens)
            st = tokens{1}{1};
            param_st = tokens{1}{2};
            param_tokens = regexp(param_st,'\s*,\s*','split');
            for i=1:numel(param_tokens)
                tk2 = regexp(param_tokens{i},'\s*(.+?)\s*=(.+)','tokens');
                phi.params.default_params.(tk2{1}{1}) = eval(tk2{1}{2});
            end
        end
        
        % parse operator
        [success, st1, st2] = parenthesisly_balanced_split(st, '<');
        
        %tokens = regexp(st, '(.+)\s*<\s*(.+)','tokens');
        if success
            phi.type='predicate';
            phi.st = st;
            phi.params.fn = [ '-(' st1 '- (' st2 '))' ];
            phi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
            return
        end
        
        [success, st1, st2] = parenthesisly_balanced_split(st, '>');
        if success
            phi.type = 'predicate';
            phi.st = st;
            phi.params.fn = [ '(' st1 ')-(' st2 ')' ];
            phi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
            return
        end
        
        % should probably discontinue this - or review its implementation
        tokens = regexp(st, '(.+)\s*=\s*(.+)', 'tokens');
        if ~isempty(tokens)
            phi.type = 'predicate';
            phi.st = st;
            if ~isfield(phi.params,'threshold')
                phi.params.threshold = 1e-14;
            end
            if ~isfield(phi.params,'max_true_value')
                phi.params.max_true_value = 1;
            end
            if ~isfield(phi.params,'alpha')
                phi.params.alpha = 1;
            end
            phi.params.fn = [ 'fun__zero(abs(' tokens{1}{1} '-(' tokens{1}{2} ')),threshold,max_true_value,alpha)'];
            phi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
            return
        end
        
        % Last possibility, the formula already exists - note: in that case
        % we ignore id and use the id of existing formula
        try
            %          id = phi.id;
            st = regexprep(st,'[()\s]','');
            phi = struct(BreachGlobOpt.STLDB(st));
            %          phi.id = id;
        catch
            error('STL_Parse',['Unknown predicate or malformed formula: ' st]);
        end
        
    case 2
        switch(varargin{1})
            case 'not'
                phi.type = 'not';
                phi.phi = varargin{2};
                
            case 'ev'
                phi.type = 'eventually';
                phi.phi = varargin{2};
                phi.interval = '[0 inf]';
                
            case 'alw'
                phi.type = 'always' ;
                phi.phi = varargin{2};
                phi.interval = '[0 inf]';
                
            case 'andn'
                phi.type = 'andn';
                phi.phin = varargin{2}; % array of STL_Formula
                
            otherwise
                phi.st = varargin{1};
                phi.evalfn = @(mode,traj,t) feval(varargin{1},mode,traj,t,varargin{2});
                phi.type = 'predicate';
        end
        
    case 3
        switch(varargin{1})
            case 'or'
                phi.type = 'or' ;
                phi.phi1 = varargin{2};
                phi.phi2 = varargin{3};
                
            case 'and'
                phi.type = 'and';
                phi.phi1 = varargin{2};
                phi.phi2 = varargin{3};
                
            case '=>'
                phi.type = '=>';
                phi.phi1 = varargin{2};
                phi.phi2 = varargin{3};
                
            case 'eventually'
                phi.type = 'eventually' ;
                phi.interval = varargin{2};
                phi.phi = varargin{3};
                
            case 'always'
                phi.type = 'always' ;
                phi.interval = varargin{2};
                phi.phi = varargin{3};
                
            case 'ev'
                phi.type = 'eventually' ;
                phi.interval = varargin{2};
                phi.phi = varargin{3};
                
            case 'alw'
                phi.type = 'always' ;
                phi.interval = varargin{2};
                phi.phi = varargin{3};
                
            case 'until'
                phi.type = 'until' ;
                phi.interval =  '[0,inf]';
                phi.phi1 = varargin{2};
                phi.phi2 = varargin{3};
                
        end
        
    case 4
        switch(varargin{1})
            case 'until'
                phi.type = 'until' ;
                phi.interval =  varargin{3};
                phi.phi1 = varargin{2};
                phi.phi2 = varargin{4};
        end
    otherwise
        error('STL_Parse','Too many arguments.')
end

end

function [success, st1, st2, interval] = parenthesisly_balanced_split(st, op)

% split st into st1 op st2 where st1 and st2 are parenthesisly balanced
success = 0;
st1 = '';
st2 = '';
interval ='';

[start_idx, end_idx, ~, ~, tokens] = regexp(st,op);

for i = 1:numel(start_idx)
    
    % checks left hand side
    st1 = st(1:start_idx(i)-1);
    st2 = st(end_idx(i)+1:end);
    
    [success, diag, st1, st2] = checks_parenthesis_balance(st1,st2);
    if success==-1
        error(['STL_Parse: expression ' st ':' diag]);
    elseif success==1
        if nargout == 4
            interval= ['[' tokens{i}{1} ']'];
        end
        st1 = strtrim(st1);
        st2 = strtrim(st2);
        
        return
    end
end

st1 = strtrim(st1);
st2 = strtrim(st2);

end

function [success, diag, st1, st2] = checks_parenthesis_balance(st1,st2)

success=0;
diag = '';

% finds parenthesis
idx_left_par1 = regexp(st1,'(');
idx_right_par1 = regexp(st1,')');

nb_left_par1 = numel(idx_left_par1);
nb_right_par1 = numel(idx_right_par1);

idx_left_par2 = regexp(st2,'(');
idx_right_par2 = regexp(st2,')');

nb_left_par2 = numel(idx_left_par2);
nb_right_par2 = numel(idx_right_par2);

% first sanity check: equal total number of ( and )

diff_par = (nb_left_par1+nb_left_par2) - (nb_right_par1+nb_right_par2);
if (diff_par>0)
    diag=sprintf('Too many (%d) opening parenthesis in expr', diff_par);
    success=-1;
    return;
elseif (diff_par<0)
    diag=sprintf('Too many (%d) closing parenthesis in expr', -diff_par);
    success=-1;
    return;
end

% checks parenthesis for st1

% first check/remove enclosing parenthesis: we should have (*par_exp where par exp is balanced
% so we check the difference in the number of left and right, if more right, then problem

diff1 = nb_left_par1- nb_right_par1; % from previous check, diff2 = -diff1

if (diff1 ~=0)
    if (diff1<0)
        success=0;
        return;
    else % alright, so diff1>0 should be the number of enclosing parenthesis
        % we remove them
        
        % checks if there is nothing but blanks before enclosing par.
        pre_st1 = st1(1:idx_left_par1(diff1));
        if ~isempty(regexp(pre_st1, '[^\(\s]'))
            success= 0;
            return;
        end
        
        % checks if there is nothing but blanks after enclosing par.
        post_st2 = st2(idx_right_par2(end-diff1+1):end);
        if ~isempty(regexp(post_st2, '[^\)\s]'))
            success= 0;
            return;
        end
        
        st1 = st1(1+idx_left_par1(diff1):end);
        idx_left_par1 = idx_left_par1(1+diff1:end);
        
        st2 = st2(1:idx_right_par2(end-diff1+1)-1);
        idx_right_par2 = idx_right_par2(1:end-diff1);
        
    end
end
% At this point, no enclosing parenthesis any more, st1 and st2 should be balanced
success = check_par(idx_left_par1, idx_right_par1) && check_par(idx_left_par2, idx_right_par2);
end

function success = check_par(idx_left_par, idx_right_par)

% idx_left_par and idx_right_par have same number of elements

assert(numel(idx_left_par) == numel(idx_right_par));
counter = 0;

%read stuff; +1 if left par, -1 if right. Whenever counter<0 exit with success==0
lcount =1;
rcount =1;
nb_par = numel(idx_left_par);

while (1)
    
    % no more left or right par.
    if (lcount > nb_par) && (rcount>nb_par)
        break;
    end
    
    if (lcount > nb_par) % no more left, add one right
        counter = counter-1;
        rcount = rcount+1;
    elseif (rcount > nb_par) % no more right, add one left
        counter=counter+1;
        lcount=lcount+1;
    else % still some left and right parenthesis
        
        next_lp = idx_left_par(lcount);
        next_rp = idx_right_par(rcount);
        if (next_lp<next_rp) % next par is left, add left
            counter=counter+1;
            lcount=lcount+1;
        else % next par is right, add right
            counter=counter-1;
            rcount=rcount+1;
        end
    end
    
    % if count went neg, no success
    if counter<0
        break;
    end
    
end

success = (counter==0); % meaning all left par have been consumed by right par

end


