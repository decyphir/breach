function NewSystem(NewSys)
%NEWSYSTEM sets up a Breach working directory for a new dynamical system
%
% Synopsis: NewSystem([NewSys])
% 
% Input:
%  The function takes zero or one argument. If no argument is provided, a
%  GUI pops up to define the new system. Otherwise, the argument should
%  be a structure which fields are described below:
%
%  - NewSys.Derivatives (cell array) must contain the dynamics as
%    strings of the form x'= x*y+x*x (expression must be a valid expression
%    in the C language)
%
%  - NewSys.x0Values (cell array) must contain the default values for
%    initial conditions of the system, as strings of the form 'x(0)=1.2'
%
%  - NewSys.ParamValues (numeric array) must contain the default values for
%    parameters of the system, as strings of the form 'param=1.2'
%
%  - NewSys.AbsTol (numeric value) if present, set the default Abs
%    tolerance
%
%  - NewSys.RelTol (numeric value) if present, set the default Rel
%    tolerance
%
%
% Example (Van der Pol system):
%   NewSys.x0Values = {'y0(0)=2.', 'y1(0)=0.'};
%   NewSys.ParamValues={'mu=2.'};
%   NewSys.Derivatives={'y0''=y1','y1''=mu*(1-y0*y0)*y1-y0'};
%   
%   NewSystem(NewSys);
%


if(nargin==0)
    
    dir = uigetdir('.','Select or create a folder for the new system');
    cd(dir);
    prompt = {'Enter derivatives (eg: x1''=a*x2*x1+b*x1*x1, x2''=x1, etc, one line per equation)', ...
        'Enter default initial conditions (x1(0)=0, x2(0)=0, etc, one line per variable)', ...
        'Enter default values for parameters (eg: a=1, b=2, etc, one line per parameter)'};
    name = 'Create new system';
    numlines = 5;
    
    default_answer = {
        sprintf(['x0''=-x1*x1-x2*x2-a*x0+a*F\n',...
                 'x1''=x0*x1-b*x0*x2-x1+G\n',...
                 'x2''=b*x0*x1+x0*x2-x2']),...
        sprintf(['x0(0)=1.0\n',...
                 'x1(0)=0.0\n',...
                 'x2(0)=0.0']),...
        sprintf(['a=0.25\n',...
                 'b=4.0\n',...
                 'F=0.5\n',...
                 'G=0.5'])...
        };
    
    answer = inputdlg(prompt,name,numlines,default_answer);
    
    if isempty(answer)
        return;
    end
    der = answer{1};
    x0val = answer{2};
    param_values = answer{3};
    
    NewSys.Derivatives={};
    NewSys.x0Values={};
    NewSys.ParamValues={};
    
    for ii = 1:size(der,1)
        new_der = strtrim(der(ii,:));
        if ~isempty(new_der)
            NewSys.Derivatives = [NewSys.Derivatives {new_der}];
        end
    end
    
    for ii = 1:size(x0val,1)
        new_x0val = strtrim(x0val(ii,:));
        if ~isempty(new_x0val)
            NewSys.x0Values = [NewSys.x0Values {new_x0val}];
        end
        
    end
    
    for ii = 1:size(param_values,1)
        new_param = strtrim(param_values(ii,:));
        if ~isempty(new_param)
            NewSys.ParamValues = [NewSys.ParamValues {new_param}];
        end
    end
    
    NewSystem(NewSys);
    
else
    x0v = NewSys.x0Values;
    deriv = NewSys.Derivatives;
    paramv = NewSys.ParamValues;
    
    variable_names = {};
    param_names = {};
    variable_values = {};
    param_values = {};
    
    % Get a name for the new system
    dr = pwd;
    indst = strfind(dr, filesep);
    Sysname = dr(indst(end)+1:end);
    
    fprintf(['\n\nSet up system ' Sysname ' in folder ' dr filesep '\n']);
    
    % writing dynamics.cpp
    
    fprintf('--- Writing file dynamics.cpp ....\n');
    fid = fopen('dynamics.cpp','w');
    
    fprintf(fid,'#include "dynamics.h"\n\n');
    
    % writing f
    
    fprintf(fid,'int f(realtype __t, N_Vector __x, N_Vector __xdot, void *__f_data) {\n\n');
    fprintf(fid,'    Fdata* __data = (Fdata*) __f_data;\n\n');
    
    dimx = numel(x0v);
    dimp = numel(paramv); %dimp is the number of system parameters
    idx_sp = [];   % indexes in paramv of the system parameters
 %  idx_pp = [];   % indexes in paramv of the properties parameter 
    %  AD: Property parameters have nothing to do in dynamics.cpp - shout if you think the contrary   

    for ii = 1:dimp
        pname = regexp(paramv{ii}, '\s*(\w+)\s*=.+', 'tokens'); % we extract the parameter name
        indexTmp = strfind(deriv,pname{1}{1});  % indexes where the current parameter appears in the equations
 %       if isempty([indexTmp{:,:}])
  %        idx_pp = [idx_pp,ii];  % if it doesn't appear, it is a property parameter
 %       else
         idx_sp = [idx_sp,ii];  % else, it is a system parameter
 %       end
    end
    dimsp = numel(idx_sp); % number of system parameters
%    dimpp = numel(idx_pp); % number of property parameters
    
    dimx_st = int2str(dimx);
    
    fprintf(fid,'    // Variables\n');
    for ii = 1:dimx % for every equation, we define the variable name
        xname = regexp(x0v{ii}, '(\w+)[ \t]*\([ \t]*(0+(\.0*)?|\.0+)[ \t]*\)\s*=\s*(.+)', 'tokens');
        variable_names = [variable_names, xname{1}{1}];
        variable_values = [variable_values, xname{1}{3}];
        fprintf(fid,['    realtype ' xname{1}{1} ' = Ith(__x,' int2str(ii-1) ');\n']);
    end
    fprintf(fid,'\n');
    
    if(dimsp>=1)
        fprintf(fid,'    // System parameters\n');
    end
    for ii = 1:dimsp % We then define the system parameters
        pname = regexp(paramv{idx_sp(ii)}, '\s*(\w+)\s*=\s*(.+)', 'tokens');
        param_names = [param_names, pname{1}{1}];
        param_values = [param_values, pname{1}{2}];
        
        fprintf(fid,['    realtype ' pname{1}{1} ' = __data->p[' int2str(ii-1+dimx) '];\n']);
    end
%    for ii = 1:dimpp % And the property parameters
%        pname = regexp(paramv{idx_pp(ii)}, '\s*(\w+)\s*=\s*(.+)', 'tokens');
%        param_names = {param_names{:}, pname{1}{1} };
%        param_values = {param_values{:}, pname{1}{2} };
%        fprintf(fid,['    realtype ' pname{1}{1} ' = data->p[' int2str(i-1+dimx+dimsp) '];\n']);
%    end
    fprintf(fid,'\n');
    
    fprintf(fid,'    // Equations\n');
    for ii = 1:dimx % we define the equations
        eqs = regexp(deriv{ii}, 'ddt{(\w+)}[ \t]*=(.+)|(\w+)''[ \t]*=(.+)', 'tokens');
         
        indexTmp = strcmp(variable_names, eqs{1}{1}); % look for the variable in x0v
        idx = find(indexTmp); %look for its index in x0v
        
        if isempty(idx)
            fclose(fid);
            error('NewSystem:equationDefinition',[' ''' eqs{1}{1} ''' has no initial value']);
        end
        if(numel(idx)>1)
            fclose(fid);
            error('NewSystem:equationDefinition',[' ''' eqs{1}{1} ''' has several initial values !?']);
        end
        
        fprintf(fid,['    Ith(__xdot,' int2str(idx-1) ')=' eqs{1}{2} ';\n' ]);
    end
    
    fprintf(fid, '\n    return(0);\n}\n');
    fprintf(fid, '\n\n');
    
    % write InitFdata and others
    
    fprintf(fid,'void InitFdata(void * &__f_data, const mxArray * mxData) {\n\n');
    fprintf(fid, '    Fdata* __data = (Fdata*) __f_data;\n\n');
    fprintf(fid,['    __data->dimx = ' dimx_st ';  // number of variables\n']);
    fprintf(fid, '    __data->dimu = 0;\n');
    fprintf(fid, '    __data->dimg = DIMG;  // number of modes\n');
    fprintf(fid,['    __data->dimp = ' int2str(dimx+dimsp) ';  //number of system parameters\n']);
    fprintf(fid, '    __data->p = mxGetPr(mxGetField(mxData,0,"p"));\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n\n');
    
    fprintf(fid, 'int UpdateFdata(realtype __t, N_Vector __x, void * &__f_data, int __Ns, N_Vector* __xS) {\n');
    fprintf(fid, '    return 0;\n}\n');
    fprintf(fid, '\n\n');
    
%    fprintf(fid, '/* gout[i] represent the distance between the system and the gard for the i-th mode */\n');
    fprintf(fid, 'int g(realtype __t, N_Vector __x, realtype *__gout, void *__f_data) {\n');
    fprintf(fid, '    return 0;\n}\n');
    fprintf(fid, '\n\n');
    
    fprintf(fid, 'int Jac(long int __N, DenseMat __J, realtype __t, \n');
    fprintf(fid, '        N_Vector __x, N_Vector __fx, void *__jac_data,\n');
    fprintf(fid, '        N_Vector __tmp1, N_Vector __tmp2, N_Vector __tmp3) {\n');
    fprintf(fid, '    return 0;\n}\n');
    
    fclose(fid);
    
    % writing dynamics.h
    
    fprintf('--- Writing file dynamics.h ....\n');
    fid = fopen('dynamics.h','w');
    
    fprintf(fid, '#ifndef DYNAMICS_H\n');
    fprintf(fid, '#define DYNAMICS_H\n\n');
    
    fprintf(fid, '#include "breach.h"\n\n');
    fprintf(fid, '#define DIMG 0 \n\n');
    
    fprintf(fid, 'class Fdata: public FdataCommon {\n');
    fprintf(fid, ' public:\n');
    fprintf(fid, '  double *p;\n');
    fprintf(fid, '};\n\n');
    
    fprintf(fid, 'int  f(realtype t, N_Vector y, N_Vector ydot, void *f_data); \n');
    fprintf(fid, 'void InitFdata(void * &f_data, const mxArray * mxData);\n');
    fprintf(fid, 'int  UpdateFdata(realtype t, N_Vector y, void * &f_data, int Ns, N_Vector* yS);\n\n');
    fprintf(fid, 'int g(realtype t, N_Vector y, realtype *gout, void *g_data);\n');
    fprintf(fid, '#endif\n');
    
    fclose(fid);
    
    % write CreateSystem.m
    
    fprintf('--- Writing file CreateSystem.m ....\n');
    fid = fopen('CreateSystem.m','w');
    
    fprintf(fid, '%% System Definition \n\n');
    fprintf(fid, 'InitBreach;\n');
    fprintf(fid, 'clear Sys;\n');
    fprintf(fid, ['Sys.DimX = ' dimx_st ';\n']);
    fprintf(fid,  'Sys.DimU = 0;\n');
    fprintf(fid, ['Sys.DimP = ' int2str(dimx+dimsp) ';\n']);
    
    fprintf(fid,'Sys.ParamList = {');
    
    for ii = 1:dimx-1
        fprintf(fid, ['''' variable_names{ii}  ''',']);
    end
    fprintf(fid, ['''' variable_names{dimx}  '''']);
    
    % The dimsp first names in param_name are the system parameters names.
    % The dimpp following ones are the property parameters names.
    for ii = 1:dimsp
        if(mod(ii,10) == 1)
            fprintf(fid, ',...\n         ');
        else
            fprintf(fid, ', ');
        end
        fprintf(fid, ['''' param_names{ii}  '''']);
    end
 %   for ii = 1:dimpp
 %       if((mod(ii,10)) == 1)
 %           fprintf(fid, ',...\n         ');
 %       else
 %           fprintf(fid, ', ');
 %       end
 %       fprintf(fid, ['''' param_names{dimsp+ii}  '''']);
 %   end
    
    fprintf(fid, '};\nSys.x0 = [');
    
    for ii = 1:dimx-1
        fprintf(fid, [variable_values{ii}  ',']);
    end
    fprintf(fid, [variable_values{dimx}  ']'';\n']);
    
    fprintf(fid, 'Sys.p = [');
    
    for ii = 1:dimx-1
        fprintf(fid, [variable_values{ii}  ',']);
    end
    fprintf(fid, [variable_values{dimx}]);
    
    for ii = 1:dimp
        if((mod(ii,10)) == 1)
            fprintf(fid, ',...\n         ');
        else
            fprintf(fid, ', ');
        end
        
        fprintf(fid,[param_values{ii}]);
    end
    fprintf(fid,']'';\n');
    
    % CVodes default options
    
    %fprintf(fid, ['Sys.Dir = ''' regexprep(dr,'\','\\\\') ''';\n']);
    fprintf(fid, 'Sys.Dir = pwd;\n');
    fprintf(fid, 'options = [];\n');
    
    if isfield(NewSys,'AbsTol')
        AbsTol = sprintf('%g,',NewSys.AbsTol);
        AbsTol = [ '[' AbsTol(1:end-1) ']' ];
        %AbsTol = num2str(NewSys.AbsTol);
    else
        AbsTol = '1.e-6';
    end
    
    if isfield(NewSys,'RelTol')
        RelTol = num2str(NewSys.RelTol);
    else
        RelTol = '1.e-6';
    end
    
    fprintf(fid, ['options = CVodeSetOptions(options,''RelTol'',' RelTol ',''AbsTol'',' AbsTol  ',''MaxNumSteps'',100000);\n']);
    fprintf(fid, 'Sys.CVodesOptions = options;\n');
    fprintf(fid, '\n%% options for sensitivity computation\n\n');
    fprintf(fid, 'pscales = ones(Sys.DimP,1);\n');
    fprintf(fid, 'FSAoptions = CVodeSetFSAOptions(''SensErrControl'',''on'',...\n');
    fprintf(fid, '                               ''ParamField'', ''p'',...\n');
    fprintf(fid, '                               ''ParamScales'',pscales );\n');
    fprintf(fid, 'Sys.CVodesSensiOptions.method = ''Simultaneous'';\n');
    fprintf(fid, 'Sys.CVodesSensiOptions.FSAoptions = FSAoptions;\n');
    
    fclose(fid);
    
    % write Init_[SysName].m
    
    fprintf(['--- Writing file Init_' Sysname '.m ....\n']);
    fid = fopen(['Init_' Sysname '.m'] ,'w');
    fprintf(fid, '%% File called previously to each trajectory computation\n');
    fprintf(fid, '%% Feel free to update it is necessary\n');
    fprintf(fid, '\n');
    fprintf(fid, ['run(''' dr filesep 'CreateSystem.m'');\n']);
    fprintf(fid, [Sysname ' = Sys;']);
    fclose(fid);
    fprintf('done.\n\n\n');
end
