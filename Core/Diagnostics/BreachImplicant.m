classdef BreachImplicant
    properties
        Intervals
        SignificantSamples
    end
    
    methods
        function obj = BreachImplicant(obj)
            obj.Intervals = [];
            obj.SignificantSamples = [];
        end 
        
        function obj = addInterval (obj, begin_value, end_value)
            if (end_value >= begin_value)
                interval.begin = begin_value;
                interval.end = end_value;
                obj.Intervals = [obj.Intervals, interval];
            end
        end
        
        function obj = addSignificantSample (obj, sample_time, sample_value)
            sample.time = sample_time;
            sample.value = sample_value;
            obj.SignificantSamples = [sample obj.SignificantSamples];
        end
        
        function intervals = getIntervals(obj)
            intervals = obj.Intervals;
        end 
        
        function interval = getInterval(obj, index)
            if (index > length(obj.Intervals))
                interval = [];
            else
                interval = obj.Intervals(index);
            end
        end 
        
        function size = getIntervalsSize(obj) 
            size = length(obj.Intervals);
        end
        
        function significant_samples = getSignificantSamples(obj)
            significant_samples = obj.SignificantSamples;
        end
    end
end