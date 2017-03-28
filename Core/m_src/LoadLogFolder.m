function B =  LoadLogFolder(folder_name)

if ~isempty(folder_name)
    dd = ls([folder_name filesep 'Br*.mat']);
    load([folder_name filesep dd(1,:)]);
    B = Br.copy();
    for id= 2:size(dd, 1)
        load([folder_name filesep dd(id,:)]);
        B.Concat(Br);
    end
    
end

end
