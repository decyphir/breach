classdef BreachTASignalGen < BreachSignalGen

    methods
        function  this = BreachTASignalGen(signals, TA_filename, labels_ranges, num_evts, time_scale)

            TA_filename = which(TA_filename);
            if isempty(TA_filename)
                error('File %s not found.', TA_filename);
            end



            % creates signal_gen and BreachSignalGen
            labels=  setdiff(fieldnames(labels_ranges), {'all__'});
            sg = TA_signal_gen2(signals,TA_filename,labels, num_evts);
            
           
            
            this = this@BreachSignalGen({sg});


            % timescale parameter
            if ~exist('time_scale', 'var')||isempty('time_scale')
                time_scale =1;
            end

            % some defaults for time steps, etc, may need adjusting
            % later
            this.SetParam('time_scale', time_scale);
            sg.min_dt = time_scale/100;
            dt = time_scale/200;
            time = 0:dt:100*time_scale;
            this.SetTime(time);

            % Adjust ranges. TODO: allow for enum and int...
            for idx_label = 1:numel(labels)
                label = labels{idx_label};


                for isig = 1:numel(signals)
                    sig = signals{isig};
                    % Assign ranges for this specific label, if any
                    if isfield(labels_ranges.(label), sig)
                        range = labels_ranges.(label).(sig);
                        params = this.expand_param_name([sig '_' label '_val']);
                        if isscalar(range) % single value
                            this.SetParam(params, range);
                        else
                            this.SetParamRanges(params, range)
                        end
                        % assign range for all labels to this label
                        % otherwise

                    elseif  isfield(labels_ranges, 'all__')&&isfield(labels_ranges.all__, sig)
                        range = labels_ranges.all__.(sig);
                        params = this.expand_param_name([sig '_' label '_val']);
                        if isscalar(range) % single value
                            this.SetParam(params, range);
                        else
                            this.SetParamRanges(params, range)
                        end
                    end


                end
            end

            % 0,1 ranges for random parameters (branching and time for
            % wordgen
            pevts = this.expand_param_name('e.*_dt');
            pbranching = this.expand_param_name('e.*_branching');

            this.SetParamRanges([pevts pbranching], [0 1]);









        end

    end

end