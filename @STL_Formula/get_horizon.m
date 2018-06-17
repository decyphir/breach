function hor = get_horizon(phi)
 % get_horizon note: uses default parameters for phi   

params = get_params(phi);
fnames = fieldnames(params);
eval_st = '';
for ip =1:numel(fnames)
    eval_st = [eval_st fnames{ip} '=' num2str(params.(fnames{ip})) ';' ];
end
eval(eval_st);

switch phi.type
    case 'predicate'
        hor=0;
    case {'always', 'eventually', 'until'}  
        interval = eval(phi.interval);
        hor_base = interval(2); 
        if hor_base == inf
            hor = inf;
            return
        end  
        children = get_children(phi );
        for ic = 1:numel(children)
            hor_children(ic) = get_horizon(children{ic});
        end
            hor = sum([hor_base hor_children]);
    otherwise
        hor_base = 0;
        children = get_children(phi );
        for ic = 1:numel(children)
            hor_children(ic) = get_horizon(children{ic});
        end
        hor = sum([hor_base hor_children]);
        
end
    

end
