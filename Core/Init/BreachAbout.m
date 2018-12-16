function BreachAbout()
% BreachAbout returns information about the breach and current system

InitBreach;
global BreachGlobOpt
breach_path = BreachGlobOpt.breach_dir;
version = BreachVersion();

disp('Breach')
disp(['Version: ' version]);
disp(['Root folder: ' breach_path]);

end