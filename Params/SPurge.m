function S = SPurge(S)
  
  % S= SPurge(S)
  %
  % Remove all fields related to a specific computation of trajectories
  %
    
  
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
  try 
    S = rmfield(S, 'props_values');
  end
  try
    S = rmfield(S, 'props');  
  end
  try
    S = rmfield(S, 'props_names');  
  end
 