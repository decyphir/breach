BrDemo.InitAFC
%%
robot = java.awt.Robot;
%%
% test to click the "Copy" button in the BreachGui
BreachGui
hBreach = findall(0, 'Tag', 'breach');
hCopyBtn = findall(hBreach, 'String', 'Copy');
posGUI = getpixelposition(hBreach);
pos = getpixelposition(hCopyBtn, true); % Get the pixel coordination relative to the GUI
jmouseemu(posGUI(1:2)+pos(1:2), 'normal')
