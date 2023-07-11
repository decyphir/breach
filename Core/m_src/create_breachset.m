function B = create_breachset(dim, varargin)
% create_breachset(dim, varargin)
%
% returns a BreachSet object of various dimensions with default names for
% parameters ('p1', 'p2', etc), with domains [0,1] 
%
%  dim is the number of parameters or a string. If it is a string, tries to
%  match it to 'load' a previous specific instance. 
% 
%  reloading can also be done using the 'load' option 
%  Bfav = create_breachset(2,'load', 'fav')

if nargin==0
    dim = 3;
end
opts.dim=dim;

if ischar(dim)
    opts.load = dim;
else
    opts.load = 'default';  % used to 'save specific previous configurations' 
end

opts.prefix_params = 'p';
opts.ranges = [0 1];
opts.samples = 1;     
opts= varargin2struct_breach(opts,varargin{:});


switch opts.load
    case 'default'
        params = arrayfun(@(c)([opts.prefix_params num2str(c)]),1:dim, 'UniformOutput',false);
        B = BreachSet(params);
        B.SetParamRanges(params,opts.ranges);
        if opts.samples > 1
            B.SampleDomain(opts.samples)
        end

end

