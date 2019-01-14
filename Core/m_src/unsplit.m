function st = unsplit(c, sep)
% unsplit convert a cell c of string into a unique string with comma
% separated elements from c

if ~exist('sep', 'var')
    sep = ',';
end

if isempty(c)
    st = '';
else
    st= c{1};
    for ic = 2:numel(c)
        st = [st sep c{ic}];
    end
end
end