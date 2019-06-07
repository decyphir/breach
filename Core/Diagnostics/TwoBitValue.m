classdef TwoBitValue
    enumeration
        FF, FT, TF, TT;
    end
    
    methods (Static)
        function [out] = getValue(v1, v2)
            if (v1 < 0 && v2 < 0)
                out = TwoBitValue.FF;
            elseif (v1 < 0 && v2 >= 0)
                out = TwoBitValue.FT;
            elseif (v1 >= 0 && v2 < 0)
                out = TwoBitValue.TF;
            else 
                out = TwoBitValue.TT;
            end                                                     
        end
    end
end