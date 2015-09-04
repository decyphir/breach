function st = short_disp(phi, maxl)
%SHORT_DISP displays a shortened version of the formula string. The firsts
% and lasts characters of the formula are keept, but the middle of the
% formula is replaced by '  \\  '.
% 
% Synopsis: st = short_disp(phi[, maxl])
% 
% Inputs:
%  - phi  : the formula we want to print
%  - maxl : (Optional, default=Inf) indicates the number of characters to
%           kepp from the beginning and the end of the formula string.
% 
% Output:
%  - st : the eventually reduced string. If get_st(phi) is shorter than
%         2*maxl+6, st is equal to get_st(phi).
%

st = disp(phi,1); %NM: why not use get_st(phi) ?

if exist('maxl','var')
    if(2*maxl+6<=length(st))
        st1 = st(1:maxl);
        st2 = st(end-maxl:end);
        st = [st1 '  \\  ' st2];
    end
end

end
