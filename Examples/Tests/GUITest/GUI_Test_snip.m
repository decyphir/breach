BrDemo.InitAFC
%%
%robot = java.awt.Robot;
%%
% test to click the "Copy" button in the BreachGui
BreachGui
hBreach = findall(0, 'Tag', 'breach');
inputEmuWrapper(hBreach, 'button_copy_set', 0.5, 'normal');
% validate
hWorkSetListBox = findall(hBreach, 'Tag', 'working_sets_listbox');
selected_name = hWorkSetListBox.String{hWorkSetListBox.Value};
added_name = hWorkSetListBox.String{end};
assert(~isempty(strfind(added_name, selected_name)))

%% test the keyboard
inputEmuWrapper(hBreach, 'edit_rename', 0.1, 'key_normal', 'testName\ENTER');

%% validate
assert(any(strcmp(hWorkSetListBox.String, 'testName')))
