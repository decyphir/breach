classdef output_gen < signal_gen

    properties % all are cells of strings, except for domains
        signals_in                  % input signals needed to compute outputs 
        domains = containers.Map()  % maps all of the above to their respective domains 
    end
    
end