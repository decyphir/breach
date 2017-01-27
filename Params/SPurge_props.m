function P = SPurge_props(P)
%SPURGE_PROPS removes all fields related to a specific computation of
% properties
%
% Synopsis:  P = SPurge_props(P)
%
% Input:
%  - P the parameter set to clean
% 
% Output:
%  - P the parameter set whithout fields related to formula evaluation.
% 
%See also SPurge SEvalProp
%

try
    P = rmfield(P, 'props_values');
end

try
    P = rmfield(P, 'props');
end

try
    P = rmfield(P, 'props_names');
end

end
