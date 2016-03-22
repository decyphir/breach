classdef BreachSet < handle
    % BreachSet Defines an API to manipulate parameters and traces. 
    %    This is the base class for BreachSystem objects. BreachSet 
    %    instances are usually also BreachSystem or some derived class 
    %    instances, so the purpose of this class is mainly to separate
    %    methods for defining and manipulate parameters and traces data.   
    % 
    % BreachSet Properties
    %   ParamRanges              - ranges of possible values for each parameter - determines the parameter sampling domain  
    %   SignalRanges             - ranges of values taken by each signal variable  
    %   AppendWhenSample=false   - when true, sampling appends new param vectors, otherwise replace.
    % 
    % BreachSet Methods
    %   GetParam            - get values for parameters, given their names 
    %   SetParam            - set values for parameters, given their names 
    %   GetParamRanges      - get ranges for parameters, given their names and ranges
    %   SetParamRanges      - set ranges for parameters, given their names and ranges
    %   GridSample          - creates parameter vectors in ParamRanges on a regularly spaced grid 
    %   CornerSample        - creates parameter vectors from the corners of ParamRanges 
    %   QuasiRandomSample   - uniformly samples n parameter vectors in ParamRanges 
    %   Concat              - concatenates parameter vectors and traces of another compatible BreachSet 
    %   PrintParams         - display the number of parameter vectors, parameter names and ranges
    %   PrintSignals        - display the number of traces, names of signal variables and ranges
    %   PlotParams          - plots parameter vectors
    %   PlotSignals         - plots signals 
    %   PlotSigPortrait     - plots signals portrait
    
    properties 
        P 
        ParamRanges  % ranges of possible values for each parameter - determines the parameter sampling domain  
        SignalRanges % ranges of values taken by each signal variable  
        AppendWhenSample=false % when true, sampling appends new param vectors, otherwise replace.
    end
        
    methods (Hidden=true)
        function P = GetP(this)
            % Get the legacy parameter set structure
            P = this.P;
        end
        
        function SetP(this, P) 
            % Get the legacy parameter set structure
            this.P = P;
        end    
    end
    
    methods
        function this = BreachSet(Sys, params, ranges)
        % BreachSet constructor from a legacy P or Sys, parameter names and ranges
            
            switch nargin
                case 0
                    return;
                case 1
                    this.P = CreateParamSet(Sys);
                case 2
                    this.P = CreateParamSet(Sys, params);
                case 3
                    this.P = CreateParamSet(Sys, params, ranges);
            end
            
            this.UpdateParamRanges();
        end
    
        function SetParam(this, params, values)
            this.P = SetParam(this.P, params, values);
            this.UpdateParamRanges();
        end
        
        function values = GetParam(this, params, ip)
            values = GetParam(this.P,params);
            if exist('ip', 'var')
                values = values(:, ip);
            end
        end
        
        % ResetParamSet Reset parameter set based on ParamRanges
        function ResetParamSet(this)
            ipr = find(diff(this.ParamRanges'));
            ranges =this.ParamRanges(ipr,:);
            if (isempty(ipr))
                this.P = CreateParamSet(this.P);
                this.P.epsi(:,:)=0;
            else
                this.P = CreateParamSet(this.P, this.P.ParamList(ipr+this.P.DimX),ranges);
            end
        end

        
        %% Get and Set param ranges
        function SetParamRanges(this, params, ranges)
            i_params = FindParam(this.P, params);
            % if we have trajectories and only set a range on property parameters, then
            % we must keep the trajectories
            save_traj = 0;
            if (isfield(this.P,'traj')&&(all(i_params>this.P.DimP)))
                save_traj=1;
                traj = this.P.traj;
                traj_ref = this.P.traj_ref;
                traj_to_compute = this.P.traj_to_compute;
                Xf = this.P.Xf;
            end
            
            this.P = SetParam(this.P,params,zeros(numel(i_params),1)); % adds new parameter name if needs be
            this.ParamRanges(i_params-this.P.DimX, :) = ranges;
            i_new_params = find(diff(this.ParamRanges'))+this.P.DimX;
            new_params = this.P.ParamList(i_new_params);
            new_ranges = this.ParamRanges(i_new_params-this.P.DimX,:);
            this.P = CreateParamSet(this.P, new_params, new_ranges);
            if (save_traj)
                this.P.traj = traj;
                this.P.traj_ref = traj_ref;
                this.P.traj_to_compute = traj_to_compute;
                this.P.Xf = Xf;
            end
        end
        
        function ranges = GetParamRanges(this, params)
            i_params = FindParam(this.P, params);
            ranges= zeros(numel(params),2);
            ranges(:,1) = -inf;
            ranges(:,2) = inf;
            for ip = 1:numel(i_params)
                if (i_params(ip)-this.P.DimX) <= size(this.ParamRanges,1)
                    ranges(ip,:) = this.ParamRanges(i_params(ip)-this.P.DimX, :);
                end
            end
        end
        
        % Set Param ranges around individual parameter vectors to zero        
        function ResetEpsi(this)
            this.P.epsi(:,:) = 0;
            this.UpdateParamRanges();
        end
        
        % Update ranges for variables from P
        function UpdateParamRanges(this)
            
            i_params = (this.P.DimX+1):numel(this.P.ParamList);
            if numel(i_params)> size(this.ParamRanges,1)
                dd = numel(i_params) - size(this.ParamRanges,1);
                mpr = size(this.ParamRanges,1);
                this.ParamRanges(mpr+1:mpr+dd,1) = inf;
                this.ParamRanges(mpr+1:mpr+dd,2) = -inf;
            end
                        
            epsis = 0*this.P.pts;
            epsis(this.P.dim,:) = this.P.epsi;
            minP = min(this.P.pts-epsis, [], 2);
            maxP = max(this.P.pts+epsis, [], 2);
            this.ParamRanges(i_params-this.Sys.DimX,:) = [minP(i_params,:), maxP(i_params,:)];
            
        end
        
        
        
        
        % Get computed trajectories
        function traces = GetTraces(this)
            traces= [];
            if isfield(this.P,'traj')
                traces = this.P.traj;
            end
        end
        
        % Update ranges for variables from trajectories in P
        function val = UpdateSignalRanges(this)
            
            if isfield(this.P, 'traj')
                if isempty(this.SignalRanges)
                    this.SignalRanges = ones(this.Sys.DimX,2);
                    minX = +inf*ones(this.Sys.DimX,1);
                    maxX = -inf*ones(this.Sys.DimX,1);
                else
                    minX = this.SignalRanges(:,1);
                    maxX = this.SignalRanges(:,2);
                end
                val=inf;
                for itraj = 1:numel(this.P.traj)
                    traj = this.P.traj(itraj);
                    traj_maxX = max(traj.X,[], 2);
                    traj_minX = min(traj.X,[], 2);
                    dist_maxX = min(maxX-traj_maxX);
                    dist_minX = min(traj_minX-minX);
                    val= min( [val dist_maxX dist_minX] );
                    minX = min([traj_minX minX],[],2);
                    maxX = max([traj_maxX maxX],[],2);
                end
                this.SignalRanges = [minX, maxX];
                this.P.SignalRanges = this.SignalRanges; % duplicate - never good I guess
            end
            
        end

        % Get signal names
        function SigNames = GetSignalNames(this)
            SigNames = this.P.ParamList(1:this.P.DimX);
        end
        
        % Plot signals
        function h = PlotSignals(this, varargin)
            if (~isfield(this.P,'traj'))
                error('No signal to plot. Use Sim command first.')
            end
      
            figure;
            h = SplotVar(this.P, varargin{:});
        end

        % Plot signals
        function h = PlotSigPortrait(this, varargin)
            if (~isfield(this.P,'traj'))
                error('No signal to plot. Use Sim command first.')
            end
      
            figure;
            SplotTraj(this.P, varargin{:});
        end
        
        % Grid Sample
        function GridSample(this, delta)
            this.P = Refine(this.P,delta);           
        end

        % Get corners of parameter domain
        function CornerSample(this)
            this.P.epsi = 2*this.P.epsi;
            this.P = Refine(this.P,2);
            this.P.epsi = this.P.epsi/2; 
        end
    
        % Quasi-Random Sample
        function QuasiRandomSample(this, delta)
            this.P = QuasiRefine(this.P,delta);
        end
    
        % Get the number of param vectors
        function nb_pts = GetNbParamVectors(this)
            nb_pts= size(this.P.pts,2);
        end
        
        % Concatenation - needs some additional compatibility checks...
        function Concat(this, other)
           this.P = SConcat(this.P, other.P); 
        end
        
        % Plot parameters
        function PlotParams(this, varargin)
            figure;
            SplotPts(this.P, varargin{:});
        end
 
        %% Printing
        function PrintSignals(this)
            if isempty(this.SignalRanges)
                disp( 'Signals:')
                disp( '-------')
                for isig = 1:this.P.DimX
                    fprintf('%s\n', this.P.ParamList{isig});
                end
            else
                
                fprintf('Signals (in range estimated over %d simulations):\n', numel(this.P.traj))
                disp('-------')
                for isig = 1:this.P.DimX
                    fprintf('%s in  [%g, %g]\n', this.P.ParamList{isig}, this.SignalRanges(isig,1),this.SignalRanges(isig,2));
                end
            end
            disp(' ')
        end

        function PrintParams(this)
            nb_pts= this.GetNbParamVectors();
            if (nb_pts==1)
                disp('Parameters:')
                disp('----------')
                for ip = this.P.DimX+1:numel(this.P.ParamList)
                    fprintf('%s=%g',this.P.ParamList{ip},this.P.pts(ip,1));
                    rg = this.ParamRanges(ip-this.P.DimX,2)-this.ParamRanges(ip-this.P.DimX,1);
                    if rg>0
                        fprintf(', can vary in [%g, %g]',this.ParamRanges(ip-this.P.DimX,1),this.ParamRanges(ip-this.Sys.DimX,2));
                    end
                    fprintf('\n');
                end
            else
                fprintf('Parameters (%d vectors):\n',nb_pts);
                disp('-------------------------');
                for ip = this.P.DimX+1:numel(this.P.ParamList)
                    rg = this.ParamRanges(ip-this.P.DimX,2)-this.ParamRanges(ip-this.P.DimX,1);
                    if rg>0
                        fprintf('%s',this.P.ParamList{ip});
                        fprintf(' varying in [%g, %g]\n',this.ParamRanges(ip-this.P.DimX,1),this.ParamRanges(ip-this.P.DimX,2));
                    else
                        fprintf('%s=%g\n',this.P.ParamList{ip},this.P.pts(ip,1));
                    end
                end
            end
            
            disp(' ')
        end
        
        %% Misc
        % Resets the system to nominal parameters
        function Reset(this)
            this.P = CreateParamSet(this.Sys);
        end
        
        % Removes computed trajectories
        function ResetSimulations(this)
            this.P = SPurge(this.P);
        end
        
        % Make a copy of a handle object - works because no property is
        % itself a handle object.
        function new = copy(this)
            % Instantiate new object of the same class.
            new = feval(class(this));
            
            % Copy all non-hidden properties.
            p = fieldnames(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end

        
    end
end