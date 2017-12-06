function [P, answer] = structdlg(P, title, prompt, default)
%STRUCTDLG is an helper function that opens a dialog window to edit the
% fields of a given structure
% 
% Synopsis: [P, answer] = structdlg(P[, title[, prompt[, default]]])
%

% manage inputs
if(nargin == 1)
    title = 'Edit params';
end
if ~exist('prompt','var')
    prompt = fieldnames(P);
end

if(~exist('default','var')||isempty(default))
    default_values = struct2cell(P);
else
    default_values = default;
end

defaultanswer = cell(1,numel(prompt));
for ii = 1:numel(prompt)
    if isnumeric(default_values{ii})
        defaultanswer{ii} = num2str(default_values{ii});
    elseif ~ischar(default_values{ii})
        defaultanswer{ii} = '';
    else
        defaultanswer{ii} = default_values{ii};
    end
end

% do stuff

options.Resize = 'on';
options.WindowStyle = 'modal';
options.Interpreter = 'None';

answer = inputdlg(prompt,title,1,defaultanswer,options);
if isempty(answer)
    P = [];
    return;
end
fnames = fieldnames(P);
for ii = 1:numel(prompt)
    %if ischar(default_values{ii})
    %  P.(fnames{i}) = answer{ii};
    %  continue;
    %end
    if(~strcmp(answer{ii},'')&&~(isempty(answer{ii})))
        try
            P.(fnames{ii}) = eval(answer{ii});
        catch %#ok<CTCH>
            P.(fnames{ii}) = answer{ii};
        end
    end
end

end
