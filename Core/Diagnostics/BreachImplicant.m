classdef BreachImplicant
    properties
        Intervals
        SignificantSamples
    end
    
    methods
        function this = BreachImplicant(this)
            this.Intervals = [];
            this.SignificantSamples = [];
        end 
        
        function this = addInterval (this, begin_value, end_value)
            if (end_value >= begin_value)
                interval.begin = begin_value;
                interval.end = end_value;
                this.Intervals = [this.Intervals, interval];
            end
        end
        
        function this = addSignificantSample (this, sample_time, sample_value)
            sample.time = sample_time;
            sample.value = sample_value;
            this.SignificantSamples = [sample this.SignificantSamples];
        end
        
        function intervals = getIntervals(this)
            intervals = this.Intervals;
        end 
        
        function interval = getInterval(this, index)
            if (index > length(this.Intervals))
                interval = [];
            else
                interval = this.Intervals(index);
            end
        end 
        
        function X = getSignal(this, time)
            X = 0*time;
            size = this.getIntervalsSize();
            for j=1:size
                interval = this.getInterval(j);
                x = interval.begin;
                y = interval.end;
                idx = (time>=x)&(time<y); 
                X(idx) = 1.;
            end
        end
        
        function size = getIntervalsSize(this) 
            size = length(this.Intervals);
        end
        
        function significant_samples = getSignificantSamples(this)
            significant_samples = this.SignificantSamples;
        end
    end
end