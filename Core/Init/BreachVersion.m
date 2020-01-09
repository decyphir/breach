function varargout = BreachVersion()
% BreachVersion returns current version of Breach

global BreachGlobOpt
VERSION_path = BreachGlobOpt.breach_dir;
v_file_name =[VERSION_path filesep 'VERSION'];
fo = fopen(v_file_name);
st = fgetl(fo);

fclose(fo);

if nargout == 0
    varargout = {};
    fprintf(st);
else
    varargout{1} = st;
end

end