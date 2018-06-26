function [ int_st ] = get_interval( phi )
%get_interval returns time interval for the top operator if any

int_st= '';
if ~(isempty(phi.interval)||isequal(phi.interval, [0,inf]))
    int_st =phi.interval;
end    

end

