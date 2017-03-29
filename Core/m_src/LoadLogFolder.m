function B =  LoadLogFolder(folder_name)
% LoadLogFolder look into folder folder_name for logged traces and returns a 
% breach system concatenating all those found. 

if ~isempty(folder_name)
    dd = dir([folder_name filesep 'Br*.mat']);
    load([folder_name filesep dd(1).name]);
    B = Br.copy();
    for id= 2:size(dd, 1)
        load([folder_name filesep dd(id).name]);
        B.Concat(Br);
    end    
end

end
