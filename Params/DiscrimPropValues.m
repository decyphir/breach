function Sets = DiscrimPropValues(S)
%
% DISCRIMPROPVALUES Separates traj set S according to properties checked for S
%   
% Synopsis: Sets = DiscrimPropValues(S)
%        
%
  
   
   if ~isfield(S,'props_values')
     Sets = S;
     return;
   end
   
   pr = S.props_values;
   nb_prop = size(pr,1);
   
   traj_status = zeros(1, size(pr,2));
   for i=1:nb_prop 
     pri = cat(1, pr(i,:).val);
     pri = (pri(:,1)>0)';
     traj_status = traj_status+pri*2.^i;
   end
          
   status = traj_status;
   
   Sets = [];
   S.status =0;
   while (~isempty(status))
     status_current = status(1);
     statusi = find(traj_status==status(1));
     status = status(status~= status(1));
     Sn = Select(S, statusi);
     Sn.status = status_current;
     Sets = [Sets Sn];
   end
   
   
function Sn = Select(S, kn)
  
  field_list= {'pts', 'epsi', 'XS0', 'Xf', 'ExpaMax', 'XSf','traj'};
  Sn = [];
  Sn.dim = S.dim;
  Sn.ParamList= S.ParamList;
  Sn.DimX = S.DimX;
  

  Sn.DimP = S.DimP;

  
  for i = 1:numel(field_list)
    
    if isfield(S, field_list{i})
      Sn.(field_list{i}) = S.(field_list{i})(:,kn);
    end
  end
  