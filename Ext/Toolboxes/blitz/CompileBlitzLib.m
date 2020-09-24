function CompileBlitzLib()

MEX = 'mex ';
INCLUDE = [ '-I.' filesep 'include '];
FLAGS = ' -c ';
src_file = ['.' filesep 'src' filesep 'libblitz.cpp '];
output_dir = ['-outdir .' filesep 'lib'];

compile_cmd = [MEX INCLUDE FLAGS src_file output_dir];

%fprintf(regexprep(compile_cmd,'\','\\\\'));
%fprintf('\n');

eval(compile_cmd);

end
