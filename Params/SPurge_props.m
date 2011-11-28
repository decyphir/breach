function S = SPurge_props(S)
% SPURGE_PROPS Removes values of satisfaction functions in a parameter set
%  
% Synopsis:  S = SPurge_props(S)
%
% Remove all fields related to a specific computation of properties
%
      
  try 
    S = rmfield(S, 'props_values');
  end

  try
    S = rmfield(S, 'props');  
  end
  
  try
    S = rmfield(S, 'props_names');  
  end
  