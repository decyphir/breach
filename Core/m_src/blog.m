function blog(st, v, v_min)
% blog(s,v,v_min) prints a string s if verbose level v is higher than some v_min

if ~exist('v_min','var')
    v_min = 0;
end
if ~exist('v','var')
    v = 0;
end

if v>=v_min
   fprintf(st); 
end

end