classdef signal_gen <handle
    % signal_gen An abstract simple class to generate single or multi-dimensional signals of some type.
    %
    % signal_gen Properties
    %    signals - cell array of names for the generated signals
    %    params  - cell array of parameter names
    %    p0      - default values for paramters, array with same dimension as params
    %
    % signal_gen Methods
    %   computeSignals (Abstract) - Takes parameter values and time and returns a signal
    %
    % See also constant_signal_gen, fixed_cp_signal_gen,
    % var_cp_signal_gen, pulse_signal_gen, step_signal_gen
    
    properties
        signals         % names of the signals it generates
        signals_domain  % domains for signals. If empty, all double, unbounded
        params          % parameters such as control points, etc
        params_domain   % domains for parameters. If empty, all double, unbounded   
        p0              % default values
    end
      
    methods
        function params = getParamNames(this)
            params = this.params;
        end
        
        function args= getSignalGenArgs(this)
            args = {};
        end
        
        function pval = get_param(this, pname)
           pidx = find(strcmp(this.params, pname));
           if ~isempty(pidx)
              pval = this.p0(pidx);  
           end
        end
        
        function set_param(this, pname, pval)
           pidx = find(strcmp(this.params, pname));
           if ~isempty(pidx)
              this.p0(pidx)=pval;  
           end
        end
    end
    
    methods (Abstract)
        computeSignals(this,p,time) % compute the signals        
    end
    
   methods
        function new = copy(this)
            % copy operator, works with R2010b or newer.
            objByteArray = getByteStreamFromArray(this);
            new = getArrayFromByteStream(objByteArray);
        end
        
        function plot(this, signal, time)
        % plot default signal  
            x =  computeSignals(this, this.p0, time);
            i_sig = find(strcmp(signal,this.signals),1);
            plot(time, x(i_sig,:), 'LineWidth', 2); 
        end
        
        function plot_enveloppe(this, signal, time, varargin)
            % Default implementation - no guarantee to generate a fair
            % enveloppe
            
            opt.max_time = 3;            
            opt.max_obj_eval = inf;
            opt = varargin2struct(opt, varargin{:});
            
            
            %%
            S = BreachSignalGen(this);                        
            dom = S.GetBoundedDomains();  
            S.SetTime(time);
            if ~isempty(dom)
                % random approach
                if 0
                    var = S.GetVariables();
                    S.CornerSample(100);
                    S.SampleDomain(var, 1000, 'quasi-random', 'append');
                    S.Sim(time);
                    hold on;
                    vcell = S.GetSignalValues(signal);
                    vbot =vcell{1};
                    vtop = vcell{2};
                    for is=  1:numel(vcell)
                        vbot = min([vbot; vcell{is}]);
                        vtop = max([vtop; vcell{is}]);
                    end
                end
                %% Optim based                
                reachmon = reach_monitor('r', signal, time);
                R = BreachRequirement(reachmon);
                pb = FalsificationProblem(S, R);
                pb.StopAtFalse = false;
                %pb.display = 'off';
                pb.log_traces = false;
                pb.max_obj_eval = opt.max_obj_eval;
                pb.max_time = opt.max_time;
                                
                pb.setup_meta();                
                pb.solve();
                [~, env] = reachmon.computeSignals(time,0,0);
                vbot = env(2,:);
                vtop = env(1,:);
                
                plot(time, vbot, 'k', 'LineWidth', .5);
                plot(time, vtop, 'k', 'LineWidth', .5);
                t2 = [time, fliplr(time)];
                inBetween = [vbot, fliplr(vtop)];
                f = fill(t2, inBetween, 'k', 'FaceAlpha', 0.1);
            end
        end                               
        
   end
    
end

