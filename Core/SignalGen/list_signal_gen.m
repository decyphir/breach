function list_sg =list_signal_gen( )

list_sg = {};

%  Get signal gen dir
InitBreach;
global BreachGlobOpt;

% gets signal_gen folder
signal_gen_dir = [ BreachGlobOpt.breach_dir filesep 'Core' filesep 'SignalGen' filesep];

% m-files
ldd = dir([signal_gen_dir '*.m']);

% each m-files might be a class
for ifile = 1:numel(ldd)
    [~, name] = fileparts(ldd(ifile).name);
    dd = superclasses(name);
    if ~isempty(dd)
        if isequal(dd{1}, 'signal_gen')
            list_sg = [list_sg name];
        end
    end
end

