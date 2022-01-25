function arg_cell = check_arg_is_cell(arg, num_args)
% check_arg_is_cell(arg, num_args); checks if arg is a cell of right size,
% or if it is a scalar. In the latter case, returns cellifyied argument. 

if nargin ==1
    if iscell(arg)
        arg_cell = arg;
        return;
    else
        num_args = 1;
    end
end

if ~iscell(arg)
        arg_cell = cell(1, num_args);
       for ia= 1:num_args 
           arg_cell{ia} = arg;
       end
else
    % already a cell checks size
    if numel(arg)~=num_args
        error('check_arg_is_cell', 'Wrong number of element in argument cell.');
    else
        arg_cell = arg;
    end
end