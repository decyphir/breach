function CompileRobusthom()
%COMPILEROBUSTNESS compiles the C++ functions for robusthom engine

global BreachGlobOpt

ext = mexext;
switch(ext)
    case {'mexw64', 'mexw32'}
        obj_ext = '.obj ';
    otherwise
        obj_ext = '.o ';
end

src_files = ['robustness.cpp ' ...
    'signal.cpp ' ...
    'mex_routines.cpp ' ...
    'av_robustness.cpp ' ...
    'window.cpp'
    ];

obj_files = ['robustness' obj_ext ...
    ' signal' obj_ext    ...
    ' mex_routines' obj_ext ...
    ' window' obj_ext ...
    ' av_robustness' obj_ext ...
    ];

MEX = 'mex ';
FLAGS = ' ';

if isfield(BreachGlobOpt, 'disable_robust_linear_interpolation')
    if BreachGlobOpt.disable_robust_linear_interpolation==1  
        FLAGS = [FLAGS '-DNO_LINEAR_INTERPOL ' ];
    end
end
    
% compiles src files
compile_cmd = [MEX FLAGS '-c  ' src_files  ];
%fprintf([regexprep(compile_cmd,'\','\\\\') '\n' ]);
eval(compile_cmd);

% compiles and robustness
compile_cmd = [MEX FLAGS '-outdir ..' filesep '.. RobustAnd.cpp ' obj_files];
%fprintf([regexprep(compile_cmd,'\','\\\\') '\n' ]);
eval(compile_cmd);

% compiles andn robustness
compile_cmd = [MEX FLAGS '-outdir ..' filesep '.. RobustAndn.cpp ' obj_files];
%fprintf([regexprep(compile_cmd,'\','\\\\') '\n' ]);
eval(compile_cmd);

% compiles or robustness
compile_cmd = [MEX FLAGS '-outdir ..' filesep '.. RobustOr.cpp ' obj_files];
%fprintf([regexprep(compile_cmd,'\','\\\\') '\n' ]);
eval(compile_cmd);

% compiles until robustness
compile_cmd = [MEX FLAGS '-outdir ..' filesep '.. RobustUntil.cpp ' obj_files];
%fprintf([regexprep(compile_cmd,'\','\\\\') '\n' ]);
eval(compile_cmd);

% compiles future robustness
compile_cmd = [MEX FLAGS '-outdir ..' filesep '.. RobustEv.cpp ' obj_files];
%fprintf([regexprep(compile_cmd,'\','\\\\') '\n' ]);
eval(compile_cmd);

% compiles average right future robustness
compile_cmd = [MEX FLAGS '-outdir ..' filesep '.. RobustAvEvRight.cpp ' obj_files];
%fprintf([regexprep(compile_cmd,'\','\\\\') '\n' ]);
eval(compile_cmd);

% compiles average left future robustness
compile_cmd = [MEX FLAGS '-outdir ..' filesep '.. RobustAvEvLeft.cpp ' obj_files];
%fprintf([regexprep(compile_cmd,'\','\\\\') '\n' ]);
eval(compile_cmd);

end
