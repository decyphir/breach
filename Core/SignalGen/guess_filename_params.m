function [params, values,pattern] = guess_filename_params(fname)
% guess_filename_params supports format
% dataname_[firstparam][val]_[secondparam][val]_[etc].[ext]

[fpath,name,ext] = fileparts(fname);

params = {};
values=[];

plist = strsplit(name,'_');

pattern = plist{1};
for ip = 2:numel(plist)
    pv = plist{ip};
    tok = regexpi(pv,'([a-z]*)(\d*)','tokens');
    if (isempty(tok)||numel(tok{1})~=2)
        error('Wrong filename format.');
    else
        params = [params tok{1}{1}];
        values = [values str2num(tok{1}{2})];
    end
    pattern= [pattern '_' tok{1}{1} '%d'];
end
pattern = [fpath '/' pattern ext];
end