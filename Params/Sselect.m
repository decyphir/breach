function Sn = Sselect(S, kn)
%SSELECT extract parameters from indices
%
% Synopsis  Pn = Sselect(P [, kn])
%
% E.g.: 
%          Peven = Sselect(P, 2:2:10)
%
%  if kn is not given, looks for a 'selected' field in S

  if (nargin==1)
    try 
      kn = find(S.selected~=0);
    catch
      Sn= S;
      return;
    end
  end
  
  Sn = [];
 
  field_list_copy = {'dim', 'ParamList', 'DimX', 'DimP', 'props_names', 'props'};

  for i = 1:numel(field_list_copy)    
    if isfield(S, field_list_copy{i})
      Sn.(field_list_copy{i}) = S.(field_list_copy{i});
    end
  end
  
  field_list_select= {'pts', 'epsi', 'selected', 'props_values'};
    
  for i = 1:numel(field_list_select)    
    if isfield(S, field_list_select{i})
      Sn.(field_list_select{i}) = S.(field_list_select{i})(:,kn);
    end
  end
  
  field_list_traj_ref_select = {'traj', 'etraj'};
  
  if isfield(S,'traj')
    if ~isfield(S, 'traj_ref')
      S.traj_ref = 1:numel();
    end
  end
 
  for i = 1:numel(field_list_traj_ref_select)    
     if isfield(S, field_list_traj_ref_select{i})
      Sn.(field_list_traj_ref_select{i}) = S.(field_list_traj_ref_select{i})(S.traj_ref(kn));
    end
  end

  field_list_traj_ref_select2 = {'XS0', 'Xf', 'ExpaMax', 'XSf'}; 
    
  for i = 1:numel(field_list_traj_ref_select2)    
     if isfield(S, field_list_traj_ref_select2{i})
      Sn.(field_list_traj_ref_select2{i}) = S.(field_list_traj_ref_select2{i})(:,S.traj_ref(kn));
    end
  end

  
 