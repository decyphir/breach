function ns = check_sensi(phis)
% checks whether properties needs apriori sensitivity computation
% TODO (for now just a quick fix for lyapunov exponent)

ns=0;

% Extract predicates

pred = [];
for ii = 1:numel(phis)
    pred = [pred STL_ExtractPredicates(phis(ii))];
    
end

% check if one of them contains lyap_exp... bof.

for ii = 1:numel(pred)
    st = disp(pred(ii), -1);
    
    if regexp(st, 'lyap_exp')
        ns = 1;
    end
end

end
