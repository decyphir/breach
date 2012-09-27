function S= SConcat(S,S2)
% SCONCAT concat two parameter sets 
% 
% Synopsis: P = SConcat(P1, P2) 
%
%  Tries to concat all compatible fields in P2 to those in P1. Basically add points and trajectories of P2 in P1. 
%
  
  
  field_list_copy = {'props_names', 'props',  'time_mult'};
 
  for i = 1:numel(field_list_copy)    
    if isfield(S2, field_list_copy{i})
      S.(field_list_copy{i}) = S2.(field_list_copy{i});
    end
  end
 
  field_list= {'pts', 'epsi', 'XS0', 'Xf', 'ExpaMax', 'XSf', ...
               'selected','props_values'};
  
  for i = 1:numel(field_list)
    
    if ((isfield(S,field_list{i}))&&(isfield(S2,field_list{i})))
      if (numel(S.(field_list{i}))==0)
        S.(field_list{i}) = S2.(field_list{i});       
      else
        S.(field_list{i}) = [S.(field_list{i}) S2.(field_list{i})];       
      end
    end    
  end

  if (isfield(S, 'traj')&&(isfield(S2,'traj')))

    %%% TEMPORARY FIX 
    if size(S.traj,1) >1
      S.traj = S.traj';
    end
    
    if size(S2.traj,1) >1
      S2.traj = S2.traj';
    end
        
    if (~isfield(S, 'traj_ref'))
      S.traj_ref= 1:numel(S.traj);
    end
    
    if (~isfield(S2, 'traj_ref'))
      S2.traj_ref= 1:numel(S2.traj);
    end
    
    S.traj_ref = [S.traj_ref S2.traj_ref+numel(S.traj)];  
    S.traj = [S.traj S2.traj];
    
  end
      
  