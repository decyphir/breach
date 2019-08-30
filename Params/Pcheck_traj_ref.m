function P = Pcheck_traj_ref(P)
% Pcheck_traj_ref(P) checks traj_ref consistency with 
        
if isfield(P,'traj')
   pts_row = P.pts(1:P.DimP,:)';
   for ir = 1:numel(P.traj_ref);
      if ~isequal(pts_row(ir,:), P.traj{P.traj_ref(ir)}.param)
          warning('inconsistent traj_ref');          
       end
   end   
end

