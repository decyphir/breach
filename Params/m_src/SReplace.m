function S = SReplace(S, Sin, kn)
  
  field_list= {'pts', 'epsi', 'XS0', 'Xf', 'ExpaMax', 'XSf','traj', ...
                      'props_values'};
    
  for i = 1:numel(field_list)    
    if isfield(S, field_list{i})
      S.(field_list{i})(:,kn) = Sin.(field_list{i});
    end
  end
