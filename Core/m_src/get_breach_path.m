function breach_path = get_breach_path() 
% get_breach_path returns breach root path  
breach_init_file = which('InitBreach');
breach_path = fileparts(breach_init_file);
end