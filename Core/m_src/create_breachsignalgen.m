function B = create_breachsignalgen(dim, varargin)
% create_breachset(dim, varargin)
%
% returns a BreachSet object of various dimensions with default names for
% parameters ('p1', 'p2', etc), with domains [0,1] 
%
%  dim is the number of parameters or a string. If it is a string, tries to
%  match it to 'load' a previous specific instance. 
% 
%  reloading can also be done using the 'load' option 
%  Bfav = create_breachsignalgen(2,'load', 'fav')

if nargin==0
    dim = 3;
end
opts.dim=dim;

if ischar(dim)
    opts.load = dim;
else
    opts.load = 'default';  % used to 'save specific previous configurations' 
end

opts.prefix_signals = 'sig';
opts.ranges = [-1 1];
opts.samples = 1;     
opts= varargin2struct_breach(opts,varargin{:});


switch opts.load
    case 'default'
        signals = arrayfun(@(c)([opts.prefix_signals num2str(c)]),1:dim, 'UniformOutput',false);
        
        sg = random_signal_gen(signals);
        B = BreachSignalGen(sg);
        
        seeds_params = B.expand_param_name('_seed');
        global_seed_param_gen = equal_param_gen('seed', seeds_params);
        B.SetParamGen(global_seed_param_gen);
            
        B.SetDomain('seed', 'int')
        if opts.samples > 1
            B.SetParam('seed',1:opts.samples)
        end

end

