function pth = BreachGetExtPath()
InitBreach;
global  BreachGlobOpt
if isfield(BreachGlobOpt, 'breach_dir')
    pth = [BreachGlobOpt.breach_dir filesep 'Ext'];
else
    pth =[];
end