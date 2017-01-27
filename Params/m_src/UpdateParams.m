function S = UpdateParams(S)
  
  field_names = fields(S.P0);
  nbf = numel(field_names);

  for i = 1:nbf
       S.pts(i) = S.P0.(field_names{i});
  end
  
 
