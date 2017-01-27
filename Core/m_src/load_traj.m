function traj = load_traj(file_name)
% LOAD_TRAJ loads a simple trace file into a traj

if nargin==0
    [FileName,PathName] = uigetfile('*.*','Load Trace File...');
    if FileName==0
        traj =[];
        return
    end
    file_name = [PathName, FileName];
end

[~, ~, ext] = fileparts(file_name);
switch (ext)
    case '.csv'
        trace = csvread(file_name,1);
    otherwise
        trace = load(file_name);
end
traj.param = trace(1,2:end);
traj.time = trace(:,1)';
traj.X = trace(:, 2:end)';

end
