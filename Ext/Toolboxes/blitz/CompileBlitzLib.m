function CompileBlitzLib()

if ismac && strcmp(computer('arch'),'maca64')
    % Workaround for recent Apple linker on macOS arm64 (maca64):
    % force classic linker for MEX to avoid undefined mex* adapter symbols.
    MEX = 'mex ''LDFLAGS=\$LDFLAGS -ld_classic'' ';
else
    MEX = 'mex ';
end
INCLUDE = [ '-I.' filesep 'include '];
FLAGS = ' -c ';
src_file = ['.' filesep 'src' filesep 'libblitz.cpp '];
output_dir = ['-outdir .' filesep 'lib'];

compile_cmd = [MEX INCLUDE FLAGS src_file output_dir];

%fprintf(regexprep(compile_cmd,'\','\\\\'));
%fprintf('\n');

eval(compile_cmd);

end
