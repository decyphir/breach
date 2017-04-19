function h_uitable = fill_uitable_params(h_uitable, params, p0, domains, issignal)
% fill_uitable_params Takes a uitable element and fill it with
% current values of parameters, domains, etc

is_signal = zeros(1, numel(p0));
if exist('issignal' , 'var')
    is_signal(issignal) = true;
end

if size(p0, 1)==1
    p0=p0';
end

set(h_uitable,'RowName',{});
set(h_uitable,'ColumnName',{'Name',  'Value', 'Range', 'Enum values', 'Type'});
set(h_uitable,'ColumnEditable', true);
set(h_uitable,'ColumnFormat', {'char', 'char', 'char', 'char', {'double','int','enum', 'bool'}});

data = cell(numel(params), 5);
for ip = 1:numel(params)
    data{ip, 1} = params{ip};
    
    if is_signal(ip)
        data{ip,2} = '(signal)';
    else
        data{ip, 2} = p0(ip);
    end
    
    if ~isempty(domains(ip).domain)
        data{ip, 3} = [ '[' num2str(domains(ip).domain) ']'] ;
    else
        data{ip, 3} = '';
    end
    data{ip, 5} = domains(ip).type;
    
    if isequal(domains(ip).type,'enum')||isequal(domains(ip).type,'int')||~isempty(domains(ip).enum)
        % try min:delta:max format
        en = domains(ip).enum;
        
        if numel(en)<=1
            data{ip, 4} = [ '  [  ' num2str(en)  '  ]'];
        else
            delta = en(2)-en(1);
            if isequal(en, en(1):delta:en(2));
                data{ip, 4} = ['  [  ' num2str(en(1)) ':' num2str(delta) ':' num2str(en(end))  '  ]'];
            else
                data{ip, 4} = ['  [  ' num2str(en) '  ]'];
            end
        end
    end
    
end
set(h_uitable, 'data', data);

end