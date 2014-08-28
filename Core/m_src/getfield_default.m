function out = getfield_default(in_struct, fieldname, default)
% GETFIELD_DEFAULT augment getfield with a default value when input structure does not have the required field
% used for handling structure options with optional fields
    
try 
    out = in_struct.(fieldname);
catch
    out = default;
end
