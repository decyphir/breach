function returnValue = fixLimit(errorCondition)
%FIXLIMIT leads the system to freeze (all derivatives to zero) if the
%errorCondition becomes true. The system must be compiled after calling
%this method. This method modifies dynamics.cpp and dynamics.h. Use this
%method only on a non-hybrid system.
%
% Synopsis : returnValue = fixLimit(errorCondition)
%
% Input :
%   - errorCondition : a string containing a valid C test. If this
%                      condition becomes true, the system is freezed.
%
% Output :
%  - returnValue : an error code. The method returns 0 if it succeeds.
%
% Example :
%
%   mySystem.Derivatives = {'ddt{x0} = k - x0'};
%   mySystem.x0Values = {'x0(0) = 1.0'};
%   mySystem.ParamValues = {'k = 2','limit = 1.8'};
%   NewSystem(mySystem);
%   fixLimit('x0<0.0 || x0>limit')
%   run(['.', filesep, 'CreateSystem']);
%   CompileSystem(Sys);
%   P = CreateParamSet(Sys,{'k'},[0 4]);
%   P = ComputeTraj(Sys,P,0:0.1:10);
%   figure;
%   SplotTraj(P);
%


% %%%%%%%%%%
% We first do the work for dynamics.cpp
% %%%%%%%%%%

fid = fopen('dynamics.cpp','r');
if(fid==-1)
    returnValue = -1;
    return ;
end

new_file = []; %used to store the modified dynamics.cpp

%we look for the function f
line = fgets(fid);
while(ischar(line) && isempty(strfind(line,' f(')))
    new_file = [new_file, line]; % we save the beginning of the file
    line = fgets(fid);
end

if ~ischar(line) % we parse the entire file and don't find the function f
    fclose(fid);
    returnValue = -2;
    return;
end

%We save the whole beginning of dynamics.cpp (including the declaration of
%the parameters and variables), untill we find the definition of the Ith
new_file = [new_file, line]; 
line = fgets(fid); %don't want to store the declaration of the function f
declaration = [];
while(ischar(line) && isempty(strfind(line,'  Ith(')))
    declaration = [declaration, line];
    new_file = [new_file, line]; % we save the beginning of the file
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -3;
    return;
end

new_file = [new_file, '    switch (__data->q) {\n    case 0 : //normal case\n'];
%**** new_file = [new_file, '    if(!(', condition, ')) {\n'];

%we copy all the definition of Ith, for both normal and degenerated cases
degenerated_case = [];
while(ischar(line) && ~isempty(strfind(line,'  Ith('))) %for each Ith
    new_file = [new_file, '    ', line]; % we write it in the new file
    degenerated_case = [degenerated_case, '    ', line(1:strfind(line,'=')), '0;\n']; %and save its definition for the degenerated case
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -4;
    return;
end

%we then add the degenerated case
new_file = [new_file, '        break;\n    case 1 : //degenerated case\n'];
new_file = [new_file, degenerated_case];
new_file = [new_file, '        break;\n    }\n'];
%**** new_file = [new_file, '    } else { //degenerated case\n'];
%**** new_file = [new_file, degenerated_case];
%**** new_file = [new_file, '    }\n'];

%Then, we look for the InitFdata function
while(ischar(line) && isempty(strfind(line,'InitFdata(')))
    new_file = [new_file, line];
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -5;
    return;
end

%we look for the last closing curly bracket
while(ischar(line) && line(1)~='}')
    new_file = [new_file, line];
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -6;
    return;
end

%just before the end of the function, we add the initialisation of the gard value
new_file = [new_file, '    __data->q = 0; // we start in the normal case\n'];

% We look for UpdateFdata
while(ischar(line) && isempty(strfind(line,'UpdateFdata(')))
    new_file = [new_file, line];
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -7;
    return;
end

new_file = [new_file, line];
new_file = [new_file, declaration];
new_file = [new_file, '    int prev_q = __data->q;\n'];
new_file = [new_file, '    if(', errorCondition, ') {\n'];
new_file = [new_file, '        __data->q = 1;\n'];
new_file = [new_file, '    }\n'];

%and we copy the end of the function until the return
line = fgets(fid);
while(ischar(line) && isempty(strfind(line,'return')))
    new_file = [new_file, line];
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -8;
    return;
end

%and we specify the return value
new_file = [new_file, '    return prev_q != __data->q;\n'];
line = fgets(fid); % we skip the previous return

%now, the function g
while(ischar(line) && isempty(strfind(line,' g(')))
    new_file = [new_file, line];
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -9;
    return;
end

% we set the value for gout[0]
new_file = [new_file, line];
new_file = [new_file, declaration];
new_file = [new_file, '    if(',errorCondition,') {\n        __gout[0] = -1;\n    } else {\n        __gout[0] = 1;\n    }\n'];

%and we copy the end of the file
line = fgets(fid);
while(ischar(line))
    new_file = [new_file, line];
    line = fgets(fid);
end

fclose(fid);
fid = fopen('dynamics.cpp','w');
fprintf(fid,new_file);
fclose(fid);

% %%%%%%%%%%
% Second step : dynamics.h
% %%%%%%%%%%

fid = fopen('dynamics.h','r');
if(fid==-1)
    returnValue = -7;
    return ;
end

% we look for the definition of DIMG
new_file = [];
line = fgets(fid);
while(ischar(line) && isempty(strfind(line, '#define DIMG')))
    new_file = [new_file, line];
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -8;
    return ;
end

% and we set DIMG at 1
new_file = [new_file, strrep(line,'0','1')]; % we replace the 0 by a 1
line = fgets(fid);

%we then look for the structure Fdata
while(ischar(line) && isempty(strfind(line, 'class Fdata')))
    new_file = [new_file, line];
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -9;
    return ;
end

%we look for the closing curly bracket
while(ischar(line) && line(1)~='}')
    new_file = [new_file, line];
    line = fgets(fid);
end

if ~ischar(line)
    fclose(fid);
    returnValue = -10;
    return;
end

%just before the end of the declaration, we add the declaration of the q
new_file = [new_file, '  int q;\n'];

%and we copy the end of the file
while(ischar(line))
    new_file = [new_file, line];
    line = fgets(fid);
end

fclose(fid);
fid = fopen('dynamics.h','w');
fprintf(fid,new_file);
fclose(fid);



returnValue = 0;

end

