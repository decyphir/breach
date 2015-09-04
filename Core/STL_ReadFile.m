function [props_names, props, signal_names] = STL_ReadFile(fname)
%STL_READFILE defines formulas from a text file.
%
% Synopsis: [props_names, prop] = STL_ReadFile(fname)
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
%          phi1 := mu until mu         # useless subformula
%          phi2 := alw_[0,2.3] phi1    # useless formula
%
%          -- end of phis.stl
%
% Outputs:
%  - props_names : is a cell array describing the names of the created
%                  formulas.
%  - props       : is a cell array containing the defined STL formulas in
%                  the same order than props_name.
%
% Example (Lorentz84):
%  [props_names,props] = STL_ReadFile('oscil_prop.stl');
%  props_names
%
%See also STL_Formula RecoverFormula
%

fid = fopen(fname,'r');

if(fid==-1)
    error('STL_ReadFile:OpeningError',['Couldn''t open file ' fname]);
end

tline = fgetl(fid);

current_id = '';
current_formula = '';
num_line = 0;
props_names = {};
props = {};
new_params = struct;

signal_names = {};
param_names = {};
p0 = [];


while ischar(tline)
    num_line = num_line+1;
    
    % first, dismiss comments  (anything starting with a #) and starting spaces
    tline = regexprep(tline, '^\s*','');
    tline = regexprep(tline, '\s*$','');
 
    if ~isempty(tline)
        % checks if we are declaring a test (from a CPSgrader spec.
        % file)
        tokens = regexp(tline, '^test\W','tokens');
        if ~isempty(tokens)
            tline='';
            break
        end
    end
   
    
    if regexp(tline, '^\#')
        tline = '';
    end
    
    tokens = regexp(tline, '\s*(.*?)\#(.+)|\s*(.*)','tokens');
    if ~isempty(tokens)
        tline = tokens{1}{1};
    else
        tline = '';
    end
    
    if ~isempty(tline)
        % checks if we are declaring signals
        tokens = regexp(tline, '^signal (.*)','tokens');
        if ~isempty(tokens)
            new_signals= strsplit(tokens{1}{1},',');
            signal_names = {signal_names{:} new_signals{:}};
            tline = '';
        end
    end
    
    if ~isempty(tline)
        % checks if we are declaring parameters
        tokens = regexp(tline, '^param (.*)','tokens');
        
        if ~isempty(current_id)
           try
                phi = STL_Formula(current_id, current_formula);
                phi = set_params(phi, new_params);
                props = [props, {phi}]; %#ok<*AGROW>              
                props_names = [props_names, {current_id}]; %#ok<AGROW>
            catch err
                fprintf(['ERROR: Problem with formula ' current_id ' at line ' ...
                   int2str(num_line-1) '\n']);
                rethrow(err);
            end
            
        end % we're done : if current_id is empty this is our first formula
        
        if ~isempty(tokens)
            param_defs = regexprep(tokens{1}{1}, '\s*','');
            param_defs = strsplit(param_defs,',');
            for ip = 1:numel(param_defs)
                name_value = strsplit(param_defs{ip},'=');
                new_params.(name_value{1})= str2num(name_value{2});
            end
            tline ='';
        end
    end
        
    if ~isempty(tline)
        % Checks if we're starting the def. of a new formula
        tokens = regexp(tline, '(\w+)\s*:=(.*)','tokens');
        if ~isempty(tokens)
            % ok try wrapping up what we have so far before starting a new formula
            if ~isempty(current_id)
                try
                    phi = STL_Formula(current_id, current_formula);
                    phi = set_params(phi, new_params);
                    props = [props, {phi}]; %#ok<*AGROW>
                    
                    props_names = [props_names, {current_id}]; %#ok<AGROW>
                catch err
                    fprintf(['ERROR: Problem with formula ' current_id ' at line ' ...
                        int2str(num_line-1) '\n']);
                    rethrow(err);
                end
                
            end % we're done : if current_id is empty this is our first formula
            
            % test the new id
            current_id = tokens{1}{1};
            try
                assignin('base', current_id, 0);
            catch %#ok<CTCH>
                error('STL_ReadFile:IdError',[current_id ' on line ' int2str(num_line) ' is not a valid id.']);
            end
            
            % start definition of formula
            current_formula = tokens{1}{2};
            
        else % we're continuing the definition of a formula
            
            if isempty(current_id)
                error('STL_ReadFile:MissingId',['On line ' int2str(num_line) ', no id given yet.']);
            end
            current_formula = [current_formula ' ' tline]; % #ok<AGROW>
        end
    end
    
    tline = fgetl(fid);
end


try
    phi = STL_Formula(current_id, current_formula);
    phi = set_params(phi, new_params);
    props = [props, {phi}]; %#ok<*AGROW>
    
    props_names = [props_names, {current_id}]; %#ok<AGROW>
catch err
    fprintf(['ERROR: Problem with formula ' current_id ' at line ' ...
        int2str(num_line-1) '\n']);
    rethrow(err);
end


end
