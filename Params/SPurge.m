function S = SPurge(S)
% SPURGE Remove all fields related to a specific computation of trajectories  
%
%  Synopsis: P = SPurge(P)
%
%  Notes: does not remove properties values
%
%  SEE ALSO SPURGE_PROPS    
  
  try
    S = rmfield(S, 'traj');
  end
  try
    S = rmfield(S, 'Xf');
  end

  try
    S = rmfield(S, 'XSf');
  end
  try
    S = rmfield(S, 'XS0');  
  end
  try
    S = rmfield(S, 'ExpaMax');  
  end  
 
 % Reset field traj_to_compute 
 
  X = S.pts(1:S.DimP,:)'; 
  [C,IA,IC] = unique(X,'rows');
  
  S.traj_ref= IC';
  S.traj_to_compute = IA';