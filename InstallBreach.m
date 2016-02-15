function InstallBreach
%INSTALLBREACH compiles C/C++ mex functions used by Breach. Needs only to
% be called once after an update of Breach

cdr = pwd;
dr = which('InstallBreach');

% setup directories of interest

breach_dir = dr(1:end-16);
breach_src_dir = [breach_dir filesep 'Core' filesep 'src']; 
qmitl_src_dir  = [breach_dir filesep '@STL_Formula' filesep 'private' filesep 'src'];

% compile qmitl C functions

MEX = 'mex ';
FLAGS = ' ';

try 
    cd(qmitl_src_dir);
    
    fprintf([MEX FLAGS '-outdir .. lemire_engine.c\n']);
    mex -outdir .. lemire_engine.c
    fprintf([MEX FLAGS '-outdir .. lemire_nd_engine.c\n']);
    mex -outdir .. lemire_nd_engine.c
    fprintf([MEX FLAGS '-outdir .. lemire_nd_maxengine.c\n']);
    mex -outdir .. lemire_nd_maxengine.c
    fprintf([MEX FLAGS '-outdir .. lemire_nd_minengine.c\n']);
    mex -outdir .. lemire_nd_minengine.c
    fprintf([MEX FLAGS '-outdir .. until_inf.c\n']);
    mex -outdir .. until_inf.c
    fprintf([MEX FLAGS '-outdir .. lim_inf.c\n']);
    mex -outdir .. lim_inf.c
    fprintf([MEX FLAGS '-outdir .. lim_inf_inv.c\n']);
    mex -outdir .. lim_inf_inv.c
    fprintf([MEX FLAGS '-outdir .. lim_inf_indx.c\n']);
    mex -outdir .. lim_inf_indx.c
    fprintf([MEX FLAGS '-outdir ../../../Core/m_src ltr.c\n']);
    mex -outdir ../../../Core/m_src/ ltr.c
    fprintf([MEX FLAGS '-outdir ../../../Core/m_src rtr.c\n']);
    mex -outdir ../../../Core/m_src/ rtr.c
    
    cd robusthom;
    CompileRobusthom;
    
catch
   error('InstallBreach:compilingSTLError','Problem compiling STL monitoring stuff.'); 
end

% compiles cvodes common stuff

sundials_dir = [breach_dir filesep 'Toolboxes' filesep 'sundials'];
sundials_inc_dir = [sundials_dir filesep 'include'];
sundials_src_dir = [sundials_dir filesep 'src' filesep 'sundials'];
sundials_cvodes_src_dir = [sundials_dir filesep 'src' filesep 'cvodes'];
sundials_nvm_src_dir = [sundials_dir filesep 'sundialsTB' filesep 'nvector' filesep 'src'];
cvodesTB_src_dir =  [breach_dir filesep 'Core' filesep 'cvodesTB++' filesep 'src'];

sundials_inc_flags = [' -I' qwrap(breach_src_dir) ...
                      ' -I' qwrap(sundials_inc_dir) ...
                      ' -I' qwrap(sundials_cvodes_src_dir) ...
                      ' -I' qwrap(cvodesTB_src_dir) ...
                      ' -I' qwrap(sundials_nvm_src_dir) ...
                    ];

sundials_src_files = [ 
           qwrap([sundials_nvm_src_dir filesep 'nvm_serial.c']) ...
           qwrap([sundials_nvm_src_dir filesep 'nvm_ops.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_band.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_bandpre.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_bbdpre.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_dense.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_diag.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodea.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_io.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodea_io.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_spils.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_spbcgs.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_spgmr.c']) ...
           qwrap([sundials_cvodes_src_dir filesep 'cvodes_sptfqmr.c']) ...
           qwrap([sundials_src_dir filesep 'sundials_band.c']) ...
           qwrap([sundials_src_dir filesep 'sundials_dense.c']) ...
           qwrap([sundials_src_dir filesep 'sundials_iterative.c']) ...
           qwrap([sundials_src_dir filesep 'sundials_nvector.c']) ...
           qwrap([sundials_src_dir filesep 'sundials_smalldense.c']) ...
           qwrap([sundials_src_dir filesep 'sundials_spbcgs.c']) ...
           qwrap([sundials_src_dir filesep 'sundials_spgmr.c']) ...
           qwrap([sundials_src_dir filesep 'sundials_sptfqmr.c']) ...
           qwrap([sundials_src_dir filesep 'sundials_math.c']) ...
           qwrap([sundials_dir filesep  'src' filesep 'nvec_ser' filesep 'nvector_serial.c' ])... 
           ];

out_dir = qwrap([breach_src_dir  filesep 'cv_obj']);

% Compose and execute compilation command for CVodes common files.
compile_cvodes = [MEX FLAGS '-c -outdir ' out_dir ' ' sundials_inc_flags ' ' sundials_src_files ];
fprintf(regexprep(compile_cvodes,'\','\\\\'));
fprintf('\n');
try 
    eval(compile_cvodes);
catch
    error('InstallBreach:compilingCVOdesError','Problem compiling CVodes');
end
% Compile blitz library

blitz_dir = [breach_dir filesep 'Toolboxes' filesep 'blitz'];
cd(blitz_dir);
try 
    CompileBlitzLib;
catch
   error('InstallBreach:compilingBlitzError','Compilation problem with Blitz++.'); 
end
% Compile mydiff

cd([breach_dir filesep 'Toolboxes' filesep 'mydiff']);
fprintf([MEX FLAGS '-lginac mydiff_mex.cpp\n']);
try
    mex -lginac mydiff_mex.cpp
catch %#ok<CTCH>
    warning('InstallBreach:noMydiffCompilation',...
            ['An error occurs when compiling mydiff. Maybe GiNaC not available.\n'...
            'Some functionnalities will not be available.']);
    fprintf('\n'); % for a nice printing
end

% cd back and clean variable
cd(cdr);
disp('-------------------------------------------------------------------');
disp('- Install successful.                                            --')
disp('- Yeah, I know: unless you really know what you''re doing, you   --');
disp('- got that "fatal error: ginac not found". This issue can be     --');
disp('- resolved under Debian-based Linux distros by installing the    --');
disp('- following packages: ginac, ginac-devel and ginac-utils.        --');
disp('- However, GinaC is used only by a few functions, so that for    --');
disp('- most users, it is still fair to say already:                   --');
disp('-------------------------------------------------------------------');
disp('- Install successful.                                            --');
disp('-------------------------------------------------------------------');

end
    
function qst = qwrap(st)
% necessary for dealing with £¨% spaces in dirnames...

qst = ['''' st ''' '];

end
