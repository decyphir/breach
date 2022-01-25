function ResetBreach(b)
% ResetBreach clear all files in ModelsData folder 

global BreachGlobOpt
d = [BreachGlobOpt.breach_dir filesep 'Ext' filesep 'ModelsData'];
if nargin<1
    choice = questdlg('Remove models data files?');
    if strcmp(choice,'Yes')
        b = 1;
    end
end

if b
  try 
      fprintf('Deleting all files and folders in %s ... \n', d);
      rmdir([d filesep 'slprj'], 's');         
      rmdir([d filesep 'ParallelTemp'], 's');
      delete([d filesep '*.*'])      
  end
end
