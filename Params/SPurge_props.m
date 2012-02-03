function S = SPurge_props(S)
% SPURGE_PROPS Removes values of satisfaction functions in a parameter set
%  
% Synopsis:  P = SPurge_props(P)
%
% Remove all fields related to a specific computation of properties
%
% SEE ALSO SPURGE
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
  