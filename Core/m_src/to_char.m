function [s, type]  = to_char(s)
% to_char(s) converts s to char from char or double

type = class(s);

if ~strcmp(type, 'char')
    try
        s = num2str(s);
    catch
        s ='';
    end
end