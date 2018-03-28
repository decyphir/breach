function compile_stl_mex(recomp, debug_level)
% This script will compile: 
%  1/ a mex function stl_eval_mex for monitoring STL formulas on arrays
%  (online monitoring) 
%  2/ a S-Function for online monitoring of STL in Simulink
%

global BreachGlobOpt

stlpp_dir = [BreachGlobOpt.breach_dir filesep 'Online'];

if ~exist('recomp', 'var')
    recomp = 0;
end

%try
%    crd = pwd;
%    cd(stlpp_dir);
%    if recomp == 1
%        system('make clean')
%    end    
%     status = system('make all');
%     if status
%        error('make not a system command.');
%     end
%     cd(crd);
%catch
    
if (~exist('debug_level','var'))
  debug_level = 0; 
end

compile_obj__=1;
compile_mex__=1;
compile_simulink__=1;
    
src_mex = 'stl_eval_mex.cpp';
src_s_function= 'onlineMonitorWrapper.cpp';

src_cpp = { ...
    'stl_parser.cpp' ...
    'stl_scanner.cpp '....
    'stl_driver.cpp '....
    'tools.cpp '....
    'stl_atom.cpp '....
    'transducer.cpp '....
    'interval_transducer.cpp '....
    'update_transducer.cpp '....
    'interval.cpp '....
    'robustness.cpp '....
    'signal.cpp '....
    'signal_expr.cpp '....
    };

%% Obj files
cxxflags = '-silent -DYYDEBUG=1 CXXFLAGS=''$CXXFLAGS -Wno-write-strings -std=gnu++11 -std=gnu++0x -Wno-deprecated-register''';
includes = ['-I' stlpp_dir filesep 'include'];
prefix_cpp = [stlpp_dir filesep 'src' filesep];
prefix_m = [stlpp_dir filesep 'src' filesep];
obj_dir = [stlpp_dir filesep 'obj' filesep];

%% bin
bin_dir = [stlpp_dir filesep 'bin' filesep];

if debug_level>0
 cxxflags = sprintf(['-DDEBUG__=%d ' cxxflags], debug_level);
end    

if ispc
 cxxflags = '-DIS_PC';
end    


if ispc
    obj_ext = 'obj';
else
    obj_ext = 'o';
end

objs= '';
for i_src = 1:numel(src_cpp)
    foo = regexprep(src_cpp(i_src),'cpp', obj_ext);
    obj = [obj_dir filesep foo{1}];
    objs = [objs ' ' obj];
end

if (compile_obj__)
    for i_src = 1:numel(src_cpp)
        cmd= sprintf('mex -c %s %s %s%s -outdir %s', cxxflags, includes,prefix_cpp, src_cpp{i_src}, obj_dir);
        if debug_level>=1
            fprintf(cmd);
            fprintf('\n');
        end
        eval(cmd);
    end
end

if (compile_mex__)
    cmd_exe = sprintf('mex  %s %s %s%s %s -outdir %s', cxxflags, includes, prefix_m, src_mex , objs, bin_dir);
    eval(cmd_exe)
end

if (compile_simulink__)
    cmd_simulink = sprintf('mex  %s %s %s%s %s -outdir %s', cxxflags, includes,prefix_m,src_s_function , objs, bin_dir);
    eval(cmd_simulink)
end
end