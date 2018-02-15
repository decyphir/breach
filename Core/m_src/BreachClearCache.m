function BreachClearCache()
% BreachClearCache clear all files in default cache location (DiskCaching)
InitBreach;
global BreachGlobOpt;
folder =[BreachGlobOpt.breach_dir filesep 'Ext' filesep 'ModelsData' filesep 'Cache'];

[status, message, messageid] = rmdir(folder, 's');
if status~=1
    if ~strcmp(messageid, 'MATLAB:RMDIR:NotADirectory')
        error(message);
    end
else
   disp(['Removed cache folder ' folder ' and all its content.']);
end
        
end