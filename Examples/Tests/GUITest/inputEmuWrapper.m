function inputEmuWrapper(handle, tag, delayI, cmd, varargin)  
%   Input parameters
%       handle - the GUI handle
%       tag - the tag name of the GUI element will be operated on 
%       delayI - delay intervers in seconds
%       cmds - cell array of the following commands{MOUSE_ACTION; KEYBOARD_ACTION 'text';}
%   Mouse action performs mouse emulation with 5 action types:
%
%      'move'|'none'|'wheel' - No mouse click (default)
%      'normal'              - Click left mouse button
%      'extend'              - Shift-click left mouse button
%      'alternate'           - Control-click left mouse button
%      'open'                - Double click any mouse button
%   INPUTEMU(KEYBOARD_ACTION,'text') performs keyboard emulation with
%   4 action types:
%
%      'key_normal'  - Normal keypress (shift-key pressed as needed)
%      'key_ctrl'    - CONTROL-key keypress
%      'key_alt'     - ALT-key keypress
%      'key_win'     - WIN-key keypress
%
%   The 'text' to be typed may contain special keys as escape characters
%   with '\' prefix:
% 
%      '\BACKSPACE'   '\TAB'         '\ENTER'       '\SHIFT'
%      '\CTRL'        '\ALT'         '\PAUSE'       '\CAPSLOCK'
%      '\ESC'         '\PAGEUP'      '\PAGEDOWN'    '\END'
%      '\HOME'        '\LEFT'        '\UP'          '\RIGHT'
%      '\DOWN'        '\PRINTSCREEN' '\INSERT'      '\DELETE'
%      '\WINDOWS'     '\NUMLOCK'     '\SCROLLLOCK'  
%
%   These escape characters are NOT case sensitive while regular characters
%   are. For regular backslash, use '\\' unless it is the only character
%   then '\' may be used. 
%
%   In addition to above action types, there are 8 low-level actions to
%   specify button or key to be down (pressed) or up (released).
%
%      'left_down'/'left_up' | 'drag_on'/'drag_off' - Left mouse button
%      'right_down'/'right_up'                      - Right mouse button
%      'middle_down'/'middle_up'                    - Middle mouse button
%      'key_down'/'key_up'                          - Keyboard key
% check the integrity of the handle and the tag

if nargin < 4
    error('Request four input variables: GUI handler and the commands');
end
if nargin == 5
    text = varargin{1};
end
    

% find unique handles & their location (its lower left hand corner)
if ~ishghandle(handle)
    error('The handle is not a figure handle');
end

hTarget = findall(handle, 'Tag', tag);
% Get the pixel coordination relative to the GUI
if isempty(hTarget)
    error(['Cannot find the Tag name' tag 'in the handle' handle]);
end
% get positions for the GUI and the tag 
posGUI = getpixelposition(handle); 
pos = getpixelposition(hTarget, true);

figure(handle);
drawnow;
switch cmd
    case {'move','none','wheel' ,'normal', 'extend', 'alternate', 'open'}
        inputemu({cmd, [posGUI(1)+pos(1),posGUI(2)+pos(2)]}')
    case {'key_normal', 'key_ctrl', 'key_alt', 'key_win'}
        inputemu({'normal', [posGUI(1)+pos(1),posGUI(2)+pos(2)]}');
        inputemu({'key_ctrl' 'a'});
        inputemu({cmd, text});
    otherwise
        error(['Command:' cmd 'does not support !'])
end
pause(delayI); % extra pause to let all mouse actions to complete

end % inputEmuWrapper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

