function l = VanderCorput(nb, base)
%  l = VanDerCorp(nb,3);
%
%  [ignore, p] = sort(l); % astuce taken from randperm.m. p is now a
%  kind of max dispersion permutation of the nb first integers


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
