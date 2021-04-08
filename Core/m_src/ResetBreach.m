function ResetBreach(b)
% ResetBreach clear all files in ModelsData folder 

global BreachGlobOpt
d = [BreachGlobOpt.breach_dir filesep 'Ext' filesep 'ModelsData'];
if nargin<1
    b=0;
end

if b||questdlg('Remove models data files?')
  try 
      rmdir([d filesep 'slprj'], 's');
      rmdir([d filesep 'ParallelTemp'], 's');
      delete([d filesep '*.*'])
  end
end
