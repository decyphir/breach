function [ st_phi, preds ] = rand_formula(n, op, preds)
%RAND_FORMULA Creates a random formula of size n, minimizing the height of the parse tree
%
% Synopsis: [ st_phi ] = rand_formula(n, op)
%
% by default op is a subset  {'and', 'ev', 'alw', 'until', 'evI', 'alwI',
% 'untilI'} (contains all operators if not specified)
%
% returns the random formula with predicates mu1, ..., mun

if nargin == 1
    op = {'or', 'and', 'ev', 'alw', 'until', 'evI', 'alwI', 'untilI'};
    preds = {'mu'};
end

if nargin == 2
    preds = {'mu'};
end

if (n==0)
    new_pred = genvarname('mu',{'mu',  preds{:}});
    preds = {'mu',preds{:}, new_pred};
    st_phi = regexprep(new_pred, 'mu(\d+)', 'x$1[t]>0' );
    
    return;
else
    opi = op{randi(numel(op))};
end

I = [0 0];
while I(1)==I(2)
    I = sort(2*rand([1 2]));
end

switch opi  
    
    case 'and'
        [st_phi1, preds] = rand_formula(floor((n-1)/2), op, preds);
        [st_phi2, preds] = rand_formula(ceil((n-1)/2), op, preds);
        st_phi  =  ['(' st_phi1  ') and (' st_phi2 ')'];
        
    case 'or'
        [st_phi1, preds] = rand_formula(floor((n-1)/2), op, preds);
        [st_phi2, preds] = rand_formula(ceil((n-1)/2), op, preds);
        st_phi  =  ['(' st_phi1  ') or (' st_phi2 ')'];       
        
    case 'not'
        [st_phi1, preds] = rand_formula(n-1, op, preds);
        st_phi  =  ['not (' st_phi1 ')'];
                
    case 'ev'
        [st_phi1, preds] = rand_formula(n-1, op, preds);
        st_phi  =  ['ev (' st_phi1 ')'];
        
    case 'alw'
        [st_phi1, preds] = rand_formula(n-1, op, preds);
        st_phi  =  ['alw (' st_phi1 ')'];
        
    case 'evI'
        [st_phi1, preds] = rand_formula(n-1, op, preds);
        st_phi  =  ['ev_[' num2str(I(1),3) ',' num2str(I(2),3)   ']  (' st_phi1 ')'];
        
    case 'once'
        [st_phi1, preds] = rand_formula(n-1, op, preds);
        st_phi  =  ['once (' st_phi1 ')'];
        
    case 'onceI'
        [st_phi1, preds] = rand_formula(n-1, op, preds);
        st_phi  =  ['once_[' num2str(I(1),3) ',' num2str(I(2),3)   ']  (' st_phi1 ')'];

    case 'hist'
        [st_phi1, preds] = rand_formula(n-1, op, preds);
        st_phi  =  ['hist (' st_phi1 ')'];
        
    case 'histI'
        [st_phi1, preds] = rand_formula(n-1, op, preds);
        st_phi  =  ['hist_[' num2str(I(1),3) ',' num2str(I(2),3)   ']  (' st_phi1 ')'];
                
    case 'alwI'
        [st_phi1, preds] = rand_formula(n-1, op, preds);
        st_phi  =  ['alw_[' num2str(I(1),3) ',' num2str(I(2),3)  ']  (' st_phi1 ')'];
        
    case 'until'
        [st_phi1, preds] = rand_formula(floor((n-1)/2), op, preds);
        [st_phi2, preds] = rand_formula(ceil((n-1)/2), op, preds);
        st_phi  =  ['(' st_phi1  ') until (' st_phi2 ')'];
        
    case 'untilI'
        [st_phi1, preds] = rand_formula(floor((n-1)/2), op, preds);
        [st_phi2, preds] = rand_formula(ceil((n-1)/2), op, preds);
        st_phi  =  ['(' st_phi1  ') until_[' num2str(I(1),3) ',' num2str(I(2),3)  '] (' st_phi2 ')'];
        
end
if (strcmp(preds{1},'mu'))
    preds = {preds{2:end}};
end
end
