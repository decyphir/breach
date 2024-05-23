function [wordgen_exe, found] = find_wordgen_exe()

Ext_path = BreachGetExtPath();
wg_path = [Ext_path filesep 'Toolboxes' filesep 'wordgen'];

if ispc
    wordgen_exe = 'wordgen.exe';
else
    wordgen_exe = 'wordgen';
end

wordgen_exe = [wg_path filesep wordgen_exe];
found = (exist(wordgen_exe, 'file')==2);

end