function  st= append_or_linebreak(st, st_new, prefix, max_char)

    if nargin <3
        prefix = '';
    end
    if nargin <4
      max_char = 100;
    end
    
      
    lines = strsplit(st,'\n');
    len_last_line = numel(lines{end});
    if len_last_line+numel(st_new)< max_char
        st = [st st_new];
    else
        st = sprintf([st '\n' prefix st_new]);
    end
        
end