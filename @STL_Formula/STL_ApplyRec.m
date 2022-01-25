function out = STL_ApplyRec(phi,func,n)
% STL_ApplyRec Recursively apply func to the phi and its subformulas

if nargin<3
    n=inf;
end

out = {func(phi)};
if n>0
    switch (phi.type)
        case 'predicate'
            return;
            
        case{'not', 'always', 'eventually', 'historically', 'once'}
            out1  = STL_ApplyRec(phi.phi,func,n-1);
            out = [out out1] ;
            
        case {'and', 'or', '=>', 'until'}
            out1  = STL_ApplyRec(phi.phi1,func,n-1);
            out2  = STL_ApplyRec(phi.phi2,func,n-1);
            out = [out out1 out2];
    end
end