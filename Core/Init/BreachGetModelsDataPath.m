function pth = BreachGetModelsDataPath()
% BreachGetExtPath Returns path to Ext folder

InitBreach;
global  BreachGlobOpt
if isfield(BreachGlobOpt, 'breach_dir')
    pth = [BreachGlobOpt.breach_dir filesep 'Ext' filesep 'ModelsData'];
else
    pth =[];
end