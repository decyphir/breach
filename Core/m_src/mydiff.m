function f_st = mydiff(f_st, param)
%MYDIFF Computes the symbolic derivative of f_st w.r.t param, where f_st is
% an expression involving symbols in param_list.
%
% Synopsis: df_st = mydiff(f_st, param)
%
% Input:
%  - f_st  : string containing the mathematical expression to differentiate
%  - param : string containing the variable w.r.t. which f_st is
%            differentiated
%
% Output:
%  - df_st : string containing the mathematical expression of
%            d f_st/d param
%
% Example:
%  
%  
%
% See also DiffFormula STL_SEvalDiff
%

f_st = regexprep(f_st,'\.\*','*');
f_st = regexprep(f_st,'\./','/');
f_st = regexprep(f_st,'\.\^','\^');

if exist('mydiff_mex','file')
    f_st = mydiff_mex(f_st,param);
    % f_st = regexprep(f_st,'\*','\.\*');
    % f_st = regexprep(f_st,'/','\./');
    % f_st = regexprep(f_st,'^','\.\^');
    return;
end

try
    syms x;
    f_st = diff(f_st,param);
catch %#ok<CTCH>
    error('mydiff:diffError','Problem: maybe Symbolic toolbox not present.');
end

%  f_st = regexprep(f_st,'\*','\.\*');
%  f_st = regexprep(f_st,'/','\./');
%  f_st = regexprep(f_st,'^','\.\^');

end
