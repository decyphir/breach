function [props_names, props, signal_names, param_names, in_signal_names, out_signal_names] = STL_ReadFile(fname, onlyLoadLastFormula)
%STL_READFILE reads formulas from a text file and loads them in the base
%workspace 
%
% Synopsis: [props_names, props, signal_names, param_names, in_signal_names, out_signal_names] = STL_ReadFile(fname)
%
% Input:
%  - fname the text file containing formulas. This text file fname should
%          consist of a sequence of definitions of the form:
%
%          phi1 := some formula
%          phi2 := some formula (might depend on phi1)
%          etc
%
%          blanks and comments beginning with a # are allowed eg :
%
%          phis.stl:
%
%          mu := x0[t] > 3             # That's a predicate
%          phi1 := mu until mu         # some subformula
%          phi2 := alw_[0,2.3] phi1    # some other formula
%
%          -- end of phis.stl
%
% Outputs:
%  - props_names : is a cell array describing the names of the created
%                  formulas.
%  - props       : is a cell array containing the defined STL formulas in
%                  the same order than props_name.
%
%See also STL_Formula 
%

% checks if Breach (and STLDB in particular) needs initialization 
global BreachGlobOpt
if isempty(BreachGlobOpt)
    evalin('base','InitBreach;');
end

fid = fopen(fname,'r');

if(fid==-1)
    error('STL_ReadFile:OpeningError',['Couldn''t open file ' fname]);
end

% JOHAN ADDED
if nargin < 2
    onlyLoadLastFormula = 0;
end
% END JOHAN ADDED

tline = fgetl(fid);

current_id = '';
current_formula = '';
num_line = 0;
props_names = {};
props = {};
new_params = struct;

signal_names = {};
in_signal_names = {};
out_signal_names = {};
param_names = {};
p0 = [];

got_it =0;

while ischar(tline)
    num_line = num_line+1;
    
    % first, dismiss comments  (anything starting with a #) and starting spaces
    tline = regexprep(tline, '^\s*','');
    tline = regexprep(tline, '\s*$','');
    
    % checks if we are declaring a test (from a CPSgrader spec. file)
    if ~isempty(tline)
        tokens = regexp(tline, '^test\W','tokens');
        if ~isempty(tokens)
            break
        end
    end
    
    % comments at start of line ignore line
    if regexp(tline, '^\#')
        tline = '';
    end
    
    % remove comment at the end of the line
    tokens = regexp(tline, '\s*(.*?)\#(.+)|\s*(.*)','tokens');
    if ~isempty(tokens)
        tline = tokens{1}{1};
    else
        tline = '';
    end
    
    
    % checks if we are declaring signals
    
    % we define here signals without I/O signatures
    if ~isempty(tline)
        tokens = regexp(tline, '^signal (.*)','tokens');
        if ~isempty(tokens)
            new_signals= strsplit(tokens{1}{1},',');
            for isig = 1:numel(new_signals)
                signal_names = {signal_names{:} strtrim(new_signals{isig})};
            end
            tline = '';
        end
    end
    
    % we define here input signals
    if ~isempty(tline)
        tokens = regexp(tline, '^input signal (.*)','tokens');
        if ~isempty(tokens)
            new_signals= strsplit(tokens{1}{1},',');
            for isig = 1:numel(new_signals)
                signal_names = {signal_names{:} strtrim(new_signals{isig})};
                in_signal_names = {in_signal_names{:} strtrim(new_signals{isig})};
            end
            tline = '';
        end
    end
    
    % we define here output signals
    if ~isempty(tline)
        tokens = regexp(tline, '^output signal (.*)','tokens');
        if ~isempty(tokens)
            new_signals= strsplit(tokens{1}{1},',');
            for isig = 1:numel(new_signals)
                signal_names = {signal_names{:} strtrim(new_signals{isig})};
                out_signal_names = {out_signal_names{:} strtrim(new_signals{isig})};
            end
            tline = '';
        end
    end
    
    % checks if we are declaring parameters
    if ~isempty(tline)
        tokens = regexp(tline, '^param (.*)','tokens');
                
        if ~isempty(tokens)
            param_defs = regexprep(tokens{1}{1}, '\s*','');
            param_defs = strsplit(param_defs,',');
            for ip = 1:numel(param_defs)
                name_value = strsplit(param_defs{ip},'=');
                new_params.(name_value{1})= str2num(name_value{2});
                param_names = unique({param_names{:} name_value{1}});
            end
            tline ='';
        end
    end
    
    % Checks if we're starting the def. of a new formula
    if ~isempty(tline)
        tokens = regexp(tline, '(\w+)\s*:=(.*)','tokens');
        if ~isempty(tokens)
            % ok try wrapping up what we have so far before starting a new formula
            if (~isempty(current_id)&& got_it == 0) && ~onlyLoadLastFormula
                try
                    phi = wrap_up(current_id, current_formula, new_params, in_signal_names, out_signal_names);
                    props = [props, {phi}]; %#ok<*AGROW>
                    props_names = [props_names, {current_id}];
                catch err
                    assignin('base', 'row_to_replace', num_line);
                    fprintf(['ERROR: Problem with formula ' current_id ' at line ' ...
                        int2str(num_line-1) '\n']);
                    rethrow(err);
                end
                
            end % we're done : if current_id is empty this is our first formula
            
            % test the new id
            current_id = tokens{1}{1};
            try
                assignin('base', current_id, 0);
                assignin('caller', current_id,0);
            catch
                error('STL_ReadFile:IdError',[current_id ' on line ' int2str(num_line) ' is not a valid id.']);
            end
            
            % start definition of formula
            current_formula = tokens{1}{2};
            got_it = 0;
        else % we're continuing the definition of a formula
            if isempty(current_id)
                error('STL_ReadFile:MissingId',['On line ' int2str(num_line) ', no id given yet.']);
            end
            current_formula = [current_formula ' ' tline]; % #ok<AGROW>
        end
    end
    
    tline = fgetl(fid);
end

fclose(fid);
try
    phi = wrap_up(current_id, current_formula, new_params, in_signal_names, out_signal_names);
    props = [props, {phi}];
    props_names = [props_names, {current_id}];
catch err   
    assignin('base', 'row_to_replace', num_line);    
    fprintf(['ERROR: Problem with formula ' current_id ' at line ' ...
        int2str(num_line-1) '\n']);
    rethrow(err);
end

end

function phi = wrap_up(current_id, current_formula, new_params, in_signal_names, out_signal_names)
phi = STL_Formula(current_id, current_formula);

% label subformulas by input/output signal names

phi = set_in_signal_names(phi, in_signal_names);
phi = set_out_signal_names(phi, out_signal_names);

% looks for potential parameters in the formula that could be overriden by
% new_params
% Note: will let function names, mex-files etc be overriden  

[~, params] = STL_ExtractSignals(phi);
% params = {};
% [~,~, ~, matches, tokens] = regexp(disp(phi,0), '(\<\w+\>)');
% for im=1:numel(matches)
%         params{end+1} = tokens{im}{1};
% end

fn = fieldnames(new_params)';
undef_params= setdiff(params, fn);
if not(isempty(undef_params)) 
    undefs = list_manip.to_string(undef_params);
    warning('STL_ReadFile:undef_params', 'Parameter(s) %s appear in formula current_id but are undeclared, default to 0. ', undefs);
end
for np = fn
   if ~any(strcmp(np{1},params))
       new_params = rmfield(new_params, np{1});
   else
       phi= set_params(phi, np{1}, new_params.(np{1}));
   end
end

assignin('base', current_id,phi);
end

