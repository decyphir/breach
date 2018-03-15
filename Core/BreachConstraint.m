classdef BreachConstraint < BreachStatus
% BreachConstraint Generic constraint class - it provides helpers to get
% data from Breach objects so different types of constraint  can easily be implemented  
    
    properties
        formula  % mostly for display (?)
        data       % data can be a BreachSet or signal or vector, etc
    end
    
    methods
        function [val, trace_vals] = Eval(this, varargin)
            % BreachConstraint.Eval returns evaluation of the constraint -
            % compute it for all traces available and returns min (implicit
            % conjunction) 
            
            % Collect traces from context
            [trajs] = this.getTraces(varargin{:});
            
            % For each trace, we collect robustness at time 0, as defined by standard semantics
            num_trajs = numel(trajs);
            trace_vals = zeros(1, num_trajs);
            for itraj = 1:num_trajs
                trace_vals(itraj) = this.evalTrace(trajs{itraj});
            end
            
            % A BreachConstraint must return a single value
            val = min(trace_vals);
            
        end
        
        
        function trajs =getTraces(this,varargin)
            %BreachSTLReq.getTraces
            
            switch numel(varargin)
                case 0   %  uses this.data
                    B = this.data;
                case 1  %  data
                    B = varargin{1};
                case 2   % time, X - easiest
                    %  orient time and X as breach usual
                    time = varargin{1};
                    X = varargin{2};
                    if size(time,2)==1
                        time=time';
                    end
                    if size(time,1) ~= 1
                        error('BreachSTLReq:data_inconsistent', 'time (first argument) should  be a one dimensional array.');
                    end
                    
                    if size(X,1)==numel(time)&&size(X,2)~=numel(time)
                        X = X';
                    end
                    if  size(X,2) ~= numel(time)
                        error('BreachSTLReq:data_inconsistent','X (second argument) should  be an array of same length as time (first arguement).');
                    end
                    
                    traj.status = 0;
                    traj.signals = this.formula.signals_in;
                    traj.params = this.formula.p0';
                    traj.time = varargin{1};
                    traj.X = varargin{2};
                    trajs = {traj};
                    if numel(traj.signals) ~= size(traj.X,1)
                        error('BreachSTLReq:data_inconsistent', 'data is inconsistent with formula' );
                    end
            end
            
            if exist('B', 'var')
                % read data , assumed to be a BreachSystem/Set
                Xs = B.GetSignalValues(this.formula.signals_in);
                if ~iscell(Xs)
                    Xs= {Xs};
                end
                
                i_params = FindParam(B.P, this.formula.params);
                
                % collect data necessary for formla evaluation
                for  i = 1:numel(Xs)
                    trajs{i}.time = B.P.traj{i}.time;
                    trajs{i}.X = Xs{i};
                    trajs{i}.signals = this.formula.signals_in;
                    try 
                        trajs{i}.params = B.P.pts(i_params, i)';   % double check this one...
                    catch
                        trajs{i}.params = this.formula.p0;
                    end
                    
                    trajs{i}.status = B.P.traj{i}.status;
                end
                
            end
        end
  
        
        
    end
    
    methods (Abstract)
        evalTrace(traj);  % specialization needs to specify what happens to traces 
     end
    
end