function expr = get_expr(phi)

if isfield(phi.params,'fn')
    expr = phi.params.fn;
else
    expr = '';
end
   
