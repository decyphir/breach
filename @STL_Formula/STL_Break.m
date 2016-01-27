function out = STL_Break(phi, n)
% STL_BREAK breaks a formula into subformulas
%
%  Synopsys:  out = STL_Break(phi, [n])
%
%  Input :
%    - phi  The formula to break.
%    - n    The depth of breaking. If not provided, the formula is broken
%             as deep as possible. If n is lower or equal to zero, an empty
%             array is returned. If n is equal to one, phi itself is
%             answered.
%
%
%  Output :
%    - out  All the subformulaes of breaking level from n to 1
%
%
%  Example :
%    STL_Formula('phi','ev ((x1[t]<1.0) and (not (x2[t]<3.0)))');
%
%    STL_Break(phi) (should) returns
%          [ x[t]<1,
%            x2[t]<3.0,
%            not (x2[t]<3.0),
%            (x1[t]<1.0) and (not (x2[t]<3.0)),
%            ev ((x1[t]<1.0) and (not (x2[t]<3.0))) ]
%
%    STL_Break(phi,2) (should) returns
%          [ (x1[t]<1.0) and (not (x2[t]<3.0)),
%            ev ((x1[t]<1.0) and (not (x2[t]<3.0))) ]
%
%    STL_Break(phi,1) (should) returns
%          [ ev ((x1[t]<1.0) and (not (x2[t]<3.0))) ]

out = [];

if nargin==1
    n = inf;
end

if n <= 0;
    return;
end

switch (phi.type)
    
    case 'predicate'
        out = phi;
        
    case {'not', 'always', 'eventually'}
        
        out = [STL_Break(phi.phi, n-1) phi] ;
        
    case {'and', 'or', '=>', 'until'}
        out = [STL_Break(phi.phi1, n-1) STL_Break(phi.phi2, n-1) phi];
end

% parameters 
params= get_params(phi);
for iphi = 1:numel(out)
   out(iphi) = set_params(out(iphi), params);
end



end

