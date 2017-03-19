% ResetBreach
global BreachGlobOpt
if questdlg('Remove models data files?')
  delete([BreachGlobOpt.breach_dir filesep 'Ext' filesep 'ModelsData' filesep '*.mat'])
  delete([BreachGlobOpt.breach_dir filesep 'Ext' filesep 'ModelsData' filesep '*.slx'])
  delete([BreachGlobOpt.breach_dir filesep 'Ext' filesep 'ModelsData' filesep '*.bck'])
end
