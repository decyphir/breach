function st = BreachVersion()
% BreachVersion returns current BreachVersion

global BreachGlobOpt
VERSION_path = BreachGlobOpt.breach_dir;
v_file_name =[VERSION_path filesep 'VERSION'];
fo = fopen(v_file_name);
st = fgetl(fo);
if nargout == 0
    disp(st);
end
fclose(fo);
end