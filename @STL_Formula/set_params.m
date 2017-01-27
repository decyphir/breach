function phi = set_params(phi, P, values)
%SET_PARAMS adds a parameter structure to a formula
%
%  syntax : phi = set_params(phi, P [, values] )
%
%  phi = set_params(phi, P) assumes P is a parameter structure of the form
%  P.p1 = val1, P.p2 = val2, etc
%
%  phi = set_params(phi, P, values) assumes P is parameter name (string)
%  or a cell of parameter names and values are corresponding values

global BreachGlobOpt

switch nargin
    
    case 2
        Pstruct =P;
    case 3
        if ischar(P)
            Pstruct = struct(P, values);
        else
            for ip = 1:numel(P)
                Pstruct.(P{ip})= values(ip);
            end
        end
         
end

fn = fieldnames(Pstruct);
if ~isempty(fn)
    for ifn= 1:numel(fn)        
        phi.params.default_params.(fn{ifn}) = Pstruct.(fn{ifn});
        % make sure the base formula gets updated with new parameters
        BreachGlobOpt.STLDB(phi.id) = phi;
    end
end
end