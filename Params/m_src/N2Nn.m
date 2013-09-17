function Nn = N2Nn(n,nb)
%N2NN Builds a set of nb^n points in N^n (nb points per axes)
%
% Synopsis : Nn = N2Nn(n,nb)
%
%  l = VanDerCorp(nb,3);
%

%  [ignore, p] = sort(l); % astuce taken from randperm.m. p is now a
% kind of max dispersion permutation of the
% nb first integers

if isscalar(nb)
    nb = nb*ones(1,n);
end
p = 1:nb(1);
Nn = N2NnIter(n,p,nb);

end

function Nn = N2NnIter(n,p,nb)
if(n==1)
    Nn = p;
else
    Nnm = N2NnIter(n-1,p,nb(1:end-1));
    nbm = size(Nnm,2);
    repmat(Nnm,1,nb(end));
    kron(1:nb(end),ones(1,nb(1)));
    Nn = [ repmat(Nnm,1,nb(end)) ; kron(1:nb(end),ones(1,nbm)) ];
end

end

function l = VanDerCorp(nb, base)

nbit = ceil(log(nb)/log(base))+1;
inv_base = 1./(base.^(1:nbit));
l = [];
for ii=0:nb-1
    l = [l baseget(ii,nbit,base)*inv_base']; %this is van der corput seq.
end

end

function b = baseget(i,nbit,base)
% i is an integer

b = zeros(1,nbit);
for kk=1:nbit
    bk = mod(i,base);
    b(kk) = bk;
    i = (i-bk)/base;
end

end
