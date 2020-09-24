classdef BreachDiagnostics  
    % BreachDiagnostics class - static utility pattern class that provides
    % functions for computing trace diagnostics
    methods (Static)
        
        function [out_implicant, error] = diag_not_f(in, in_implicant, samples)
            %DIAG_NOT_F explains why a not phi formula is false
            %
            %  syntax : [out_implicant, error] = diag_not_f(in, in_implicant, samples)
            %
            %  INPUT:
            %  in - robustness signal for phi in the form in.times in.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why not phi is false
            %  error - 0 if there was no error encountered, 1 otherwise
            [out_implicant, error] = BreachDiagnostics.diag_unary_plogic (in, in_implicant, samples);
        end  
        
        function [out_implicant, error] = diag_not_t(in, in_implicant, samples)
            %DIAG_NOT_T explains why a not phi formula is true
            %
            %  syntax : [out_implicant, error] = diag_not_t(in, in_implicant, samples)
            %
            %  INPUT:
            %  in - robustness signal for phi in the form in.times in.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why not phi is true
            %  error - 0 if there was no error encountered, 1 otherwise
            [out_implicant, error] = BreachDiagnostics.diag_unary_plogic (in, in_implicant, samples);
        end    
        
        function [out1_implicant, out2_implicant, error] = diag_or_f(in1, in2, in_implicant, samples)
            %DIAG_OR_F explains why a phi1 or phi2 formula is false
            %
            %  syntax : [out1_implicant, out2_implicant, error] = diag_or_f(in1, in2, in_implicant, samples)
            %
            %  INPUT:
            %  in1 - robustness signal for phi1 in the form in1.times in1.values
            %  in2 - robustness signal for phi2 in the form in2.times in2.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why phi1 or phi2 is false
            %  error - 0 if there was no error encountered, 1 otherwise
            
            %  phi1 or phi2 is false throughout the interval [a,b]
            %  because both phi1 and phi2 are false throughout [a,b]
            out1_implicant = BreachImplicant;
            out2_implicant = BreachImplicant;
            
            % Copy the intervals
            intervals = in_implicant.getIntervals();
            for(i=1:length(intervals))
                interval = intervals(i);
                out1_implicant = out1_implicant.addInterval(interval.begin, interval.end);
                out2_implicant = out1_implicant.addInterval(interval.begin, interval.end);
            end
            
            % Update the significant samples
            samples = in_implicant.getSignificantSamples();
            for(i=1:length(samples))
                sample = samples(i);
                v1 = interp1(in1.times, in1.values, sample.time, 'previous');
                v2 = interp1(in2.times, in2.values, sample.time, 'previous');
                
                out1_implicant.addSignificantSample(sample.time, v1);
                out2_implicant.addSignificantSample(sample.time, v2);
            end
            
            error = 0;
        end
        
        function [out1_implicant, out2_implicant, error] = diag_or_t(in1, in2, in_implicant, samples)
            %DIAG_OR_T explains why a phi1 or phi2 formula is true
            %
            %  syntax : [out1_implicant, out2_implicant, error] = diag_or_f(in1, in2, in_implicant, samples)
            %
            %  INPUT:
            %  in1 - robustness signal for phi1 in the form in1.times in1.values
            %  in2 - robustness signal for phi2 in the form in2.times in2.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why phi1 or phi2 is true
            %  error - 0 if there was no error encountered, 1 otherwise
            error = 0;
            
            % phi1 or phi2 is true throughout the interval [a,b] because
            % either phi1 is true only
            % or phi2 is true only
            % or both phi1 and phi2 are true
            % 
            % we partition [a,b] into intervals that each covers an above
            % case, and update the explanations accordingly
            % The actual explanation is done in diag_binary_plogic
            [out1_implicant, out2_implicant] = ...
                BreachDiagnostics.diag_binary_plogic(BreachOperator.OR, in1, in2, in_implicant, samples);
        end
        
        function [out1_implicant, out2_implicant, error] = diag_and_f(in1, in2, in_implicant, samples)
            %DIAG_AND_F explains why a phi1 and phi2 formula is false
            %
            %  syntax : [out1_implicant, out2_implicant, error] = diag_and_f(in1, in2, in_implicant, samples)
            %
            %  INPUT:
            %  in1 - robustness signal for phi1 in the form in1.times in1.values
            %  in2 - robustness signal for phi2 in the form in2.times in2.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why phi1 and phi2 is false
            %  error - 0 if there was no error encountered, 1 otherwise
            
            % phi1 and phi2 is false throughout the interval [a,b] because
            % either phi1 is false only
            % or phi2 is false only
            % or both phi1 and phi2 are false
            % 
            % we partition [a,b] into intervals that each covers an above
            % case, and update the explanations accordingly
            % The actual explanation is done in diag_binary_plogic
            
            error = 0;
            [out1_implicant, out2_implicant] = ...
                BreachDiagnostics.diag_binary_plogic(BreachOperator.AND, in1, in2, in_implicant, samples);
        end    
        
        function [out1_implicant, out2_implicant, error] = diag_and_t(in1, in2, in_implicant, samples)
            %DIAG_AND_T explains why a phi1 and phi2 formula is true
            %
            %  syntax : [out1_implicant, out2_implicant, error] = diag_and_t(in1, in2, in_implicant, samples)
            %
            %  INPUT:
            %  in1 - robustness signal for phi1 in the form in1.times in1.values
            %  in2 - robustness signal for phi2 in the form in2.times in2.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why phi1 and phi2 is true
            %  error - 0 if there was no error encountered, 1 otherwise
            
            %  phi1 and phi2 is false throughout the interval [a,b]
            %  because both phi1 and phi2 are true throughout [a,b]
            out1_implicant = BreachImplicant;
            out2_implicant = BreachImplicant;
            
            % Copy the intervals
            intervals = in_implicant.getIntervals();
            for(i=1:length(intervals))
                interval = intervals(i);
                out1_implicant = out1_implicant.addInterval(interval.begin, interval.end);
                out2_implicant = out1_implicant.addInterval(interval.begin, interval.end);
            end
            
            % Update the significant samples
            samples = in_implicant.getSignificantSamples();
            for(i=1:length(samples))
                sample = samples(i);
                v1 = interp1(in1.times, in1.values, sample.time, 'previous');
                v2 = interp1(in2.times, in2.values, sample.time, 'previous');
                
                out1_implicant.addSignificantSample(sample.time, v1);
                out2_implicant.addSignificantSample(sample.time, v2);
            end
            
            error = 0;
        end
             
        function [out1_implicant, out2_implicant, error] = diag_implies_f(in1, in2, in_implicant, samples)
            %DIAG_IMPLIES_F explains why phi1 => phi2 formula is false
            %
            %  syntax : [out1_implicant, out2_implicant, error] = diag_implies_f(in1, in2, in_implicant, samples)
            %
            %  INPUT:
            %  in1 - robustness signal for phi1 in the form in1.times in1.values
            %  in2 - robustness signal for phi2 in the form in2.times in2.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why phi1 => phi2 is false
            %  error - 0 if there was no error encountered, 1 otherwise
            
            %  phi1 => phi2 is false throughout the interval [a,b]
            %  because both phi1 is true and phi2 is false throughout [a,b]
            out1_implicant = BreachImplicant;
            out2_implicant = BreachImplicant;
            
            % Copy the intervals
            intervals = in_implicant.getIntervals();
            for(i=1:length(intervals))
                interval = intervals(i);
                out1_implicant = out1_implicant.addInterval(interval.begin, interval.end);
                out2_implicant = out2_implicant.addInterval(interval.begin, interval.end);
            end
            
            % Update the significant samples
            samples = in_implicant.getSignificantSamples();
            for(i=1:length(samples))
                sample = samples(i);
                v1 = interp1(in1.times, in1.values, sample.time, 'previous');
                v2 = interp1(in2.times, in2.values, sample.time, 'previous');
                
                if(sample.value == -v1)
                    out1_implicant = out1_implicant.addSignificantSample(sample.time, v1);
                else
                    out2_implicant = out2_implicant.addSignificantSample(sample.time, v2);
                end
            end
            
            error = 0;
        end
        
        function [out1_implicant, out2_implicant, error] = diag_implies_t(in1, in2, in_implicant, samples)
            %DIAG_IMPLIES_T explains why a phi1 and phi2 formula is true
            %
            %  syntax : [out1_implicant, out2_implicant, error] = diag_implies_t(in1, in2, in_implicant, samples)
            %
            %  INPUT:
            %  in1 - robustness signal for phi1 in the form in1.times in1.values
            %  in2 - robustness signal for phi2 in the form in2.times in2.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why phi1 => phi2 is true
            %  error - 0 if there was no error encountered, 1 otherwise
            
            % phi1 => phi2 is true throughout the interval [a,b] because
            % either phi1 is false and phi2 is true
            % or phi1 is false (while phi2 is false)
            % or phi2 is true (while phi2 is true)
            % 
            % we partition [a,b] into intervals that each covers an above
            % case, and update the explanations accordingly
            % The actual explanation is done in diag_binary_plogic
            error = 0;
            [out1_implicant, out2_implicant] = ...
                BreachDiagnostics.diag_binary_plogic(BreachOperator.IMPLIES, in1, in2, in_implicant, samples);
        end
      
        function [out_implicant, error] = diag_alw_f(in, bound, in_implicant, samples)
            %DIAG_ALW_F explains why a alw_I phi is false
            %
            %  syntax : [out_implicant, error] = diag_alw_f(in, bound, in_implicant, samples)
            %
            %  INPUT:
            %  in - robustness signal for phi in the form in.times in.values
            %  bound - interval I that bounds the always operator
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why alw_I phi is false
            %  error - 0 if there was no error encountered, 1 otherwise
            
            %  alw_I phi is false throughout the intervals [a,b] because
            %  of all the segments where phi is false in the interval
            %  [a,b]
            %  These actual segments are computed by diag_unary_tlogic
            [out_implicant, error] = BreachDiagnostics.diag_unary_tlogic(BreachOperator.ALW, in, bound, in_implicant, samples);
        end
        
        function [out_implicant, error] = diag_alw_t(in, bound, in_implicant, samples)
            %DIAG_ALW_T explains why alw_I phi is true
            %
            %  syntax : [out_implicant, error] = diag_alw_t(in, bound, in_implicant, samples)
            %
            %  INPUT:
            %  in - robustness signal for phi in the form in.times in.values
            %  bound - interval I that bounds the always operator
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why alw_I phi is true
            %  error - 0 if there was no error encountered, 1 otherwise
            
            %  alw_I phi is true throughout the intervals [a,b] because
            %  phi is true throughout [a,b] + I
            out_implicant = BreachImplicant;
            error = 0;
            size = in_implicant.getIntervalsSize();
            for(i=1:size)
                interval = in_implicant.getInterval(i);
                new_begin = interval.begin + bound.begin;
                new_end = interval.end + bound.end;
                out_implicant = out_implicant.addInterval(new_begin, new_end);
                
                interval.begin = new_begin;
                interval.end = new_end;
                out_tmp = BreachDiagnostics.diag_signal_restrict_to_interval(in, interval);
                for(j=1:length(samples))
                    for(k=1:length(out_tmp.values))
                        if(out_tmp.values(k) == samples(j).value)
                            out_implicant = out_implicant.addSignificantSample(out_tmp.times(k), out_tmp.values(k));
                        end
                    end
                end
            end
            
            
        end
        
        function [out_implicant, error] = diag_ev_f(in, bound, in_implicant, samples)
            %DIAG_EV_F explains why ev_I phi is false
            %
            %  syntax : [out_implicant, error] = diag_ev_f(in, bound, in_implicant, samples)
            %
            %  INPUT:
            %  in - robustness signal for phi in the form in.times in.values
            %  bound - interval I that bounds the always operator
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why ev_I phi is false
            %  error - 0 if there was no error encountered, 1 otherwise
            
            %  ev_I phi is true throughout the intervals [a,b] because
            %  phi is false throughout [a,b] + I
            out_implicant = BreachImplicant;
            error = 0;
            size = in_implicant.getIntervalsSize();
            for(i=1:size)
                interval = in_implicant.getInterval(i);
                new_begin = interval.begin + bound.begin;
                new_end = interval.end + bound.end;
                out_implicant = out_implicant.addInterval(new_begin, new_end);
                
                interval.begin = new_begin;
                interval.end = new_end;
                out_tmp = BreachDiagnostics.diag_signal_restrict_to_interval(in, interval);
                for(j=1:length(samples))
                    for(k=1:length(out_tmp.values))
                        if(out_tmp.values(k) == samples(j).value)
                            out_implicant = out_implicant.addSignificantSample(out_tmp.times(k), out_tmp.values(k));
                        end
                    end
                end
       
            end
        end
        
        function [out_implicant, error] = diag_ev_t(in, bound, in_implicant, samples)
            %DIAG_EV_T explains why a ev_I phi is true
            %
            %  syntax : [out_implicant, error] = diag_ev_t(in, bound, in_implicant, samples)
            %
            %  INPUT:
            %  in - robustness signal for phi in the form in.times in.values
            %  bound - interval I that bounds the always operator
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why ev_I phi is true
            %  error - 0 if there was no error encountered, 1 otherwise
            
            %  ev_I phi is true throughout the intervals [a,b] because
            %  of all the segments where phi is true in the interval
            %  [a,b]
            %  These actual segments are computed by diag_unary_tlogic
            [out_implicant, error] = BreachDiagnostics.diag_unary_tlogic(BreachOperator.EV, in, bound, in_implicant, samples);
        end
    end
    
    methods(Static,Access=private)
        function [out_implicant, error] = diag_unary_plogic(in, in_implicant, samples)
            %DIAG_UNARY_PLOGIC private function for explaining why not phi
            % is true or false
            %
            %  syntax : [out_implicant, error] = diag_unary_plogic(in, in_implicant, samples)
            %
            %  INPUT:
            %  in - robustness signal for phi in the form in.times in.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why not phi is true or false
            %  error - 0 if there was no error encountered, 1 otherwise
            
            %  not phi is true (false) throughout the intervals [a,b] because
            %  phi is (false) true throughout the interval [a,b]
            out_implicant = BreachImplicant;
            
            % Copy the intervals
            intervals = in_implicant.getIntervals();
            for(i=1:length(intervals))
                interval = intervals(i);
                out_implicant = out_implicant.addInterval(interval.begin, interval.end);
            end
            
            % Update the significant samples
            samples = in_implicant.getSignificantSamples();
            for(i=1:length(samples))
                sample = samples(i);
                
                out_implicant = out_implicant.addSignificantSample(sample.time, -sample.value);
            end
            
            error = 0;
        end
        
        function [out_implicant, error] = diag_unary_tlogic(operator, in, bound, in_implicant, samples)
            %DIAG_UNARY_TLOGIC private function that explains why ev_I phi/alw_I phi is
            %true/false
            %
            %  syntax : [out_implicant, error] = diag_unary_tlogic(in, bound, in_implicant, samples)
            %
            %  INPUT:
            %  in - robustness signal for phi in the form in.times in.values
            %  bound - interval I that bounds the always operator
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why ev_I phi/alw_I is true/false
            %  error - 0 if there was no error encountered, 1 otherwise
            
            %  ev_I phi is true throughout the intervals [a,b] because
            %  phi is false throughout [a,b] + I
            out_implicant = BreachImplicant;
            error = 0;
            size = in_implicant.getIntervalsSize();
            % We want to explain every interval in the in_implicant
            for(i=1:size)
                interval_to_explain = in_implicant.getInterval(i);
                flag = 1;
                % After each iteration, we explain part of
                % interval_to_explain, hence the interval_to_explain
                % becomes smaller after each iteration. We compute 
                % explanations until we have explained the entire
                % interval
                while (interval_to_explain.end >= interval_to_explain.begin && flag)
                    if (interval_to_explain.end == interval_to_explain.begin)
                        flag = 0;
                    end
                    % We want to explain an interval [a,b] for a formula
                    % ev_[c,d] phi/alw_[c,d] phi. Hence we should search
                    % for explanations in the interval [a+c,b+d]
                    search_interval.begin = interval_to_explain.begin + bound.begin;
                    search_interval.end = interval_to_explain.end + bound.end;
                    % We restrict the input signal 
                    % to the segment [a+c, b+d]
                    tmp = BreachDiagnostics.diag_signal_restrict_to_interval(in, search_interval);
                    tmp_size = length(tmp.times);
                    % We fetch the first sample
                    prev_time = tmp.times(1);
                    prev_value = tmp.values(1);
                    begin_time = inf;
                    end_time = inf;
                    % The first sample satisfies phi and the operator is
                    % ev_I phi, or it violates phi and the operator is
                    % alw_I phi -- we have found a beginning of an
                    % explanation
                    if (prev_value >= 0 && operator == BreachOperator.EV)
                        begin_time = prev_time;
                    elseif (prev_value < 0 && operator == BreachOperator.ALW)
                        begin_time = prev_time;    
                    end
                
                    for j=2:tmp_size
                        current_time = tmp.times(j);
                        current_value = tmp.values(j);
                        
                        % The current sample satisfies phi, the previous sample violates phi
                        % and the operator is ev_I phi, or
                        % the current sample violates phi, the previous sample satisfies phi 
                        % and the operator is alw_I phi -- 
                        % we have found a beginning of an explanation
                        if (current_value >= 0 && prev_value < 0 && operator == BreachOperator.EV)
                            begin_time = current_time;
                        elseif (current_value < 0 && prev_value >= 0 && operator == BreachOperator.ALW)
                            begin_time = current_time;
                        % The current sample violates phi, the previous sample satisfies phi
                        % and the operator is ev_I phi, or
                        % the current sample satisfies phi, the operator is
                        % ev_I and we are at the end of the search window,
                        % or
                        % the current sample satisfies phi, the previous sample violates phi 
                        % and the operator is alw_I phi, or
                        % the current sample violates phi, the operator is
                        % alw_I and we are at the end of the search window,
                        % -- we have found the end of an explanation
                        elseif (((current_value < 0 && prev_value >= 0) || ...
                                (current_value >= 0 && j == tmp_size)) && operator == BreachOperator.EV)
                            end_time = current_time;
                            out_implicant = out_implicant.addInterval(begin_time, end_time);
                        elseif (((current_value >= 0 && prev_value < 0) || ...
                                (current_value < 0 && j == tmp_size)) && operator == BreachOperator.ALW)
                            end_time = current_time;
                            out_implicant = out_implicant.addInterval(begin_time, end_time);
                        end
                        samples = in_implicant.getSignificantSamples();
                        for k=1:length(samples)
                            sample = samples(k);
                            if (sample.value == prev_value)
                                out_implicant = out_implicant.addSignificantSample(prev_time, prev_value);
                            end
                        end
                        
                        prev_time = current_time;
                        prev_value = current_value;
                    end
                    
                    
                    % Find how much we have explained, and remove it from
                    % interval_to_explain
                    interval_to_explain.begin = end_time - bound.begin;
                end
                error = 0;
            end
        end
        
        function [out1_implicant, out2_implicant] = diag_binary_plogic(operator, in1, in2, in_implicant, samples)
            %DIAG_BINARY_PLOGIC private function that explains why (1) phi1
            %and phi2 formula is false, (2) phi1 or phi2 is true or (3)
            %phi1 => phi2 is true
            %
            %  syntax : [out1_implicant, out2_implicant, error] = diag_binary_plogic(in1, in2, in_implicant, samples)
            %
            %  INPUT:
            %  in1 - robustness signal for phi1 in the form in1.times in1.values
            %  in2 - robustness signal for phi2 in the form in2.times in2.values
            %  in_implicant - BreachImplicant containing intervals to be explained 
            %  samples - robustness samples that guides the explanation
            %  OUTPUT:
            %  out_implicant - BreachImplicant containing intervals
            %  explaining why (1) phi1 and phi2 is false, (2) phi1 or phi2 is true, or (3) phi1 => phi2 is true
            %  error - 0 if there was no error encountered, 1 otherwise
            
            out1_implicant = BreachImplicant;
            out2_implicant = BreachImplicant;
            
            size = in_implicant.getIntervalsSize();
            % We want an explanation for each interval in the in_implicant
            for i = 1:size
                interval = in_implicant.getInterval(i);
                
                % We first restrict the search in signals in1 and in2
                % to the interval [a,b] to explain
                out1_tmp = BreachDiagnostics.diag_signal_restrict_to_interval(in1, interval);
                out2_tmp = BreachDiagnostics.diag_signal_restrict_to_interval(in2, interval);
                % We normalze the two signals so that they have the same 
                % samples
                [out1_tmp, out2_tmp] = BreachDiagnostics.diag_two_signals_normalize_sampling (out1_tmp, out2_tmp);
                
                % We pick the first sample from restricted input signals
                % and build a two-bit value:
                % FF - if v1 < 0 and v2 < 0
                % FT - if v1 < 0 and v2 >= 0
                % TF - if v1 >= 0 and v2 < 0
                % TT - if v1 >= 0 and v2 >= 0
                begin_time = out1_tmp.times(1);
                begin_idx = 1;
                old_v1 = out1_tmp.values(1);
                old_v2 = out2_tmp.values(1);
                old_value = TwoBitValue.getValue(old_v1,old_v2);
                
                if length(out1_tmp.times)==1
                    end_time = out1_tmp.times(1);
                    end_idx = 1;
                    new_v1 = out1_tmp.values(1);
                    new_v2 = out2_tmp.values(1);
                    
                    new_value = TwoBitValue.getValue(new_v1,new_v2);
                    [out1_implicant, out2_implicant] = BreachDiagnostics.diag_binary_plogic_update_implicants( ...
                            operator, old_value, out1_implicant, out2_implicant, begin_time, end_time, begin_idx, end_idx, out1_tmp, out2_tmp, samples);
                    
                else
                    for j=2:length(out1_tmp.times)
                        t = out1_tmp.times(j);
                        new_v1 = out1_tmp.values(j);
                        new_v2 = out2_tmp.values(j);
                        
                        new_value = TwoBitValue.getValue(new_v1,new_v2);
                        
                        if (new_value ~= old_value || j == length(out1_tmp.times))
                            end_time = t;
                            end_idx = j;
                            [out1_implicant, out2_implicant] = BreachDiagnostics.diag_binary_plogic_update_implicants( ...
                                operator, old_value, out1_implicant, out2_implicant, begin_time, end_time, begin_idx, end_idx, out1_tmp, out2_tmp, samples);
                            begin_time = t;
                            begin_idx = j;
                            old_value = new_value;
                            old_v1 = new_v1;
                            old_v2 = new_v2;
                        end
                    end
                end
            end
        end   
        
        function [out1, out2] = diag_binary_plogic_update_implicants(oper, value, in1, in2, btime, etime, bidx, eidx, s1, s2, samples)
            out1 = in1;
            out2 = in2;
            
            switch(oper)
                case BreachOperator.OR
                    if (value == TwoBitValue.FT)
                        out2 = out2.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    elseif (value == TwoBitValue.TF)
                        out1 = out1.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    elseif (value == TwoBitValue.TT)
                        flag1 = 0;
                        flag2 = 0;
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                    flag1 = 1;
                                end
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                    flag2 = 1;
                                end
                            end
                        end
                        if (flag1)
                           out1 = out1.addInterval(btime, etime);
                        end
                        if (flag2)
                           out2 = out2.addInterval(btime, etime);
                        end
                    end                
                case BreachOperator.AND
                    if (value == TwoBitValue.FT)
                        out1 = in1.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    elseif (value == TwoBitValue.FF)
                        flag1 = 0;
                        flag2 = 0;
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                    flag1 = 1;
                                end
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                    flag2 = 1;
                                end
                            end
                        end
                        if (flag1)
                           out1 = out1.addInterval(btime, etime);
                        end
                        if (flag2)
                           out2 = out2.addInterval(btime, etime);
                        end
                    elseif (value == TwoBitValue.TF)
                        out2 = in2.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    end
                case BreachOperator.IMPLIES
                    if (value == TwoBitValue.TT)
                        out2 = in2.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    elseif (value == TwoBitValue.FT)
                        flag1 = 0;
                        flag2 = 0;
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == -sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                    flag1 = 1;
                                end
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                    flag2 = 1;
                                end
                            end
                        end
                        if (flag1)
                           out1 = out1.addInterval(btime, etime);
                        end
                        if (flag2)
                           out2 = out2.addInterval(btime, etime);
                        end
                    elseif (value == TwoBitValue.FF)
                        out1 = in1.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == -sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    end
            end
        end
        
        function [out1, out2] = diag_binary_plogic_update_implicants_with_both_io(oper, value, in1, in2, btime, etime, bidx, eidx, s1, s2, samples)
            out1 = in1;
            out2 = in2;
            
            switch(oper)
                case BreachOperator.OR
                    if (value == TwoBitValue.FT)
                        out2 = out2.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    elseif (value == TwoBitValue.TF)
                        out1 = out1.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    elseif (value == TwoBitValue.TT)
                        flag1 = 0;
                        flag2 = 0;
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                    flag1 = 1;
                                end
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                    flag2 = 1;
                                end
                            end
                        end
                        if (flag1)
                           out1 = out1.addInterval(btime, etime);
                        end
                        if (flag2)
                           out2 = out2.addInterval(btime, etime);
                        end
                    end
                
                case BreachOperator.AND
                    if (value == TwoBitValue.FT)
                        out1 = in1.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    elseif (value == TwoBitValue.FF)
                        flag1 = 0;
                        flag2 = 0;
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                    flag1 = 1;
                                end
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                    flag2 = 1;
                                end
                            end
                        end
                        if (flag1)
                           out1 = out1.addInterval(btime, etime);
                        end
                        if (flag2)
                           out2 = out2.addInterval(btime, etime);
                        end
                    elseif (value == TwoBitValue.TF)
                        out2 = in2.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    end
                case BreachOperator.IMPLIES
                    if (value == TwoBitValue.TT)
                        out2 = in2.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    elseif (value == TwoBitValue.FT)
                        flag1 = 0;
                        flag2 = 0;
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == -sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                    flag1 = 1;
                                end
                                if (s2.values(i) == sample.value)
                                    out2 = out2.addSignificantSample(sample.time, sample.value);
                                    flag2 = 1;
                                end
                            end
                        end
                        if (flag1)
                           out1 = out1.addInterval(btime, etime);
                        end
                        if (flag2)
                           out2 = out2.addInterval(btime, etime);
                        end
                    elseif (value == TwoBitValue.FF)
                        out1 = in1.addInterval(btime, etime);
                        for(i=bidx:eidx)
                            for(j=1:length(samples))
                                sample = samples(j);
                                if (s1.values(i) == -sample.value)
                                    out1 = out1.addSignificantSample(sample.time, sample.value);
                                end
                            end
                        end
                    end
            end
        end
        
        function [out] = diag_signal_restrict_to_interval(in, interval)
            t = in.times;
            v = in.values;
            itv = [interval.begin, interval.end];
            
            t_tmp = union(t, itv,'sorted');
            v_tmp = interp1(t, v, t_tmp, 'previous');
            
            size = length(t_tmp);
            t_out = [];
            v_out = [];
            for (i=1:size)
                if(t_tmp(i) >= interval.begin && t_tmp(i) <= interval.end)
                    t_out = [t_out t_tmp(i)];
                    v_out = [v_out v_tmp(i)];
                end
            end
            
            out.times = t_out;
            out.values = v_out;     
        end    
        
        function [in1, in2] = diag_two_signals_normalize_sampling(in1, in2)
            t1 = in1.times;
            v1 = in1.values;
            t2 = in2.times;
            v2 = in2.values;
            
            t_tmp = union(t1, t2, 'sorted');
            
            if ~isequal(t_tmp,t1)
                v1_tmp = interp1(t1, v1, t_tmp, 'previous');
                in1.times = t_tmp;
                in1.values = v1_tmp;
            end
            
            if ~isequal(t_tmp,t2)
                v2_tmp = interp1(t2, v2, t_tmp, 'previous');
                in2.times = t_tmp;
                in2.values = v2_tmp;
            end
            
        end

    end
end