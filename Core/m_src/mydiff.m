function st = mydiff(st, param)
% function df_st = mydiff(f_st, param, param_list)
%
%  Computes the symbolic derivative of f_st w.r.t p, where f_st is an
%  expression involving symbols in param_list.
%

st = regexprep(st,'\.\*','*');
st = regexprep(st,'\./','/');
st = regexprep(st,'\.\^','\^');

if exist('mydiff_mex') % NM: replace by exist(('mydiff_mex','class') ?
    st = mydiff_mex(st,param);
    % st = regexprep(st,'\*','\.\*');
    % st = regexprep(st,'/','\./');
    % st = regexprep(st,'^','\.\^');
    return;
end

try
    syms x;
    st = diff(st,param);
catch
    error('mydiff:diffError','Problem: maybe Symbolic toolbox not present.');
end

%  st = regexprep(st,'\*','\.\*');
%  st = regexprep(st,'/','\./');
%  st = regexprep(st,'^','\.\^');

end
