function dphi = DiffFormula(phi, param)
%DIFFFORMULA Symbolically derives all predicates of a formula with respect
% to a parameter
%
% Synopsis: dphi = DiffFormula(phi, param)
%
% Input:
%  - phi   : the formula to derivate. phi must not contain ddt{} or d{}{}
%            structure, or an error is thrown.
%  - param : one parameter name
%
% Output:
%  - dphi : a STL_Formula. The name of this STL_Structure is 'd'
%           followed by the name of phi
%
% See also STL_SEvalDiff
%

switch (phi.type)
    
    case 'predicate'
        dphi = DiffPredicate(phi,param);
        
    case 'not'
        dphi1 = DiffFormula(phi.phi,param);
        dphi = STL_Formula(['d' phi.id],'not', dphi1);
        
    case 'or'
        dphi1 = DiffFormula(phi.phi1,param);
        dphi2 = DiffFormula(phi.phi2,param);
        dphi = STL_Formula(['d' phi.id], 'or',dphi1, dphi2);
        
    case 'and'
        dphi1 = DiffFormula(phi.phi1,param);
        dphi2 = DiffFormula(phi.phi2,param);
        dphi = STL_Formula(['d' phi.id], 'and',dphi1, dphi2);
        
    case '=>'
        dphi1 = DiffFormula(phi.phi1,param);
        dphi2 = DiffFormula(phi.phi2,param);
        dphi = STL_Formula(['d' phi.id], '=>',dphi1,  dphi2);
        
    case 'always'
        dphi1 = DiffFormula(phi.phi,param);
        dphi = STL_Formula(['d' phi.id], 'always', phi.interval, dphi1);
        
    case 'eventually'
        dphi1 = DiffFormula(phi.phi,param);
        dphi = STL_Formula(['d' phi.id], 'eventually', phi.interval, dphi1);
        
    case 'until'
        dphi1 = DiffFormula(phi.phi1,param);
        dphi2 = DiffFormula(phi.phi2,param);
        dphi = STL_Formula(['d' phi.id], 'until', dphi1, phi.interval, dphi2);
end

end

function dphi = DiffPredicate(PHI, param)

dphi = PHI;
dphi.id = ['d' PHI.id];

% parse param predicates

st = PHI.st;

% on regarde s'il y a des paramètre de propriété à l'ancienne mode
% (défini dans la formule par " | param=1" )
tokens = regexp(st, '(.+)\s*\|\s*(.+)','tokens');

if ~isempty(tokens)
    st = tokens{1}{1};
    param_st = tokens{1}{2};
    param_tokens = regexp(param_st,'\s*,\s*','split');
    for i=1:numel(param_tokens)
        tk2 = regexp(param_tokens{i},'\s*(.+?)\s*=(.+)','tokens');
        dphi.params.(tk2{1}{1}) = eval(tk2{1}{2});
    end
end

% parse operator

tokens = regexp(st, '(.+)\s*<\s*(.+)','tokens');
if ~isempty(tokens)
    
    dphi.type='predicate';
    
    lhs = tokens{1}{1};
    lhs = DiffMu(lhs, param);
    rhs = tokens{1}{2};
    rhs = DiffMu(rhs, param);
    
    dphi.st = [lhs ' < ' rhs];
    dphi.params.fn = [ '-(' lhs  '-' rhs ')' ]; % NM: why not [rhs '-' lhs] AD:good question, but there must be a reason
    dphi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
    return
    
end

tokens = regexp(st, '(.+)\s*>\s*(.+)','tokens');
if ~isempty(tokens)
    dphi.type='predicate';
    
    lhs = tokens{1}{1};
    lhs = DiffMu(lhs, param);
    rhs = tokens{1}{2};
    rhs = DiffMu(rhs, param);
    
    dphi.st = [lhs ' > ' rhs];
    
    dphi.params.fn = [ lhs '-' rhs ];  
    dphi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
    return
end

tokens = regexp(st, '(.+)\s*=\s*(.+)', 'tokens');
if ~isempty(tokens)
    dphi.type='predicate';
    
    
    lhs = tokens{1}{1};
    lhs = DiffMu(lhs, param);
    
    rhs = tokens{1}{2};
    rhs = DiffMu(rhs, param);
    dphi.st = [lhs ' = ' rhs];
    
    if ~isfield(PHI.params,'threshold')
        dphi.params.threshold = 1e-14;
    end
    if ~isfield(PHI.params,'max_true_value')
        dphi.params.max_true_value = 1;
    end
    if ~isfield(PHI.params,'alpha')
        dphi.params.alpha = 1;
    end
    dphi.params.fn = [ 'fun__zero(abs(' lhs '-' rhs '),threshold,max_true_value,alpha)'];
    dphi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
    return
end

end

function dmu = DiffMu(mu, param)
%DIFFMU returns a string containing the analyticial expression of
% dmu/dparam

if strcmp(mu,'0')
    dmu = mu;
else
    
    fn = mu;
    
    % checks for double differentiation
    
    [~,~,~,matches] = regexp(fn, ['ddt\{\s*' '\w+' '\s*\}\[(.+?)\]']);
    if (numel(matches)>0)
        error('DiffFormula:doubleDifferentiation',...
                ['\nThe formula already contains derivatives of variables, ',...
                'cannot differentiate more than once.\n',...
                'You can resolve this by replacing "ddt{your_var}[t]" by ',...
                '"(your_var[t+epsi]-your_var[t])/epsi" where epsi is ',...
                'small wrt the evolution time scale of your_var.\n']);
    end
    
    [~,~,~,matches] = regexp(fn, ['d\{\s*' '\w+' '\s*\}{(.+?)}\[(.+?)\]']);
    if (numel(matches)>0)
        error('DiffFormula:doubleDifferentiation',...
                'The formula already contains derivatives of variables, cannot differentiate more than once')
    end
    
    % find occurrences of variables that we will need to diff against
    
    [~,~,~,matches,tokens] = regexp(fn, '(\w+?)\[(.+?)\]');
    
    nm = numel(matches);
    vartmp = cell(1,nm);
    variables = cell(1,nm);
    times = cell(1,nm);
    for jj = 1:nm
        variables{jj} = tokens{jj}{1};
        times{jj} = tokens{jj}{2};
        
        % replaces xj[tj] with xj__tj_ in the expression
        
        vartmp{jj} = [variables{jj} '__t' num2str(jj) '_'];
        fn = regexprep(fn, regexptranslate('escape', [variables{jj} '[' times{jj} ']']), vartmp{jj},'once');
        
    end
    
    % From now on, we implement the formula
    %
    %   dfn(x1,..,xn,p)/dp  = dfn/dp + dfn/dx1 dx1/dp + .. + dfn/dxn dxn/dp
    %                          \                                        /
    %                            -----------------  -------------------
    %                                             \/
    %
    %                                     partial derivatives
    
    
    %  call to mydiff for param and each vartmp{j}
    
    dfn = mydiff(fn, param); % dfn/dp (partial)
    
    if strcmp(dfn,'0')
        dfn = '';
    end
    
    for jj=1:nm
        var = vartmp{jj};
        dfn_dxj = mydiff(fn, var);
        
        dxj_dp = ['d{' variables{jj} '}{' param '}__t' num2str(jj) '_' ];
        
        % sum results, multiplied by dxj_dp
        
        if strcmp(dfn_dxj,'1')
            if isempty(dfn)
                dfn = dxj_dp;
            else
                dfn = [dfn '+' dxj_dp]; %#ok<AGROW>
            end
        elseif ~strcmp(dfn_dxj,'0')
            if isempty(dfn)
                dfn = ['(' dfn_dxj ')*' dxj_dp];
            else
                dfn = [dfn '+(' dfn_dxj ')*' dxj_dp]; %#ok<AGROW>
            end
        end
    end
    
    if strcmp(dfn,'')
        dmu = '0';
        return;
    end
    
    % replace time specifications
    for jj=1:nm
        dfn = regexprep(dfn,['__t' num2str(jj) '_'],[ '[' times{jj} ']']);
    end
    dmu = dfn;
    
    %  TODO ? - simplify (?) -> mysimplify.cpp using ginac
    
end

end
