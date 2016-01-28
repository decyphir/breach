classdef BreachSet < handle
    properties
        P % legacy Breach structure, contains data and trajectories, etc
        ParamRanges
        SignalRanges
        AppendWhenSample=false 
    end
    
    methods
        
        % Get and Set of the parameter set structure
        function P = GetP(this)
            P = this.P;
        end
        function SetP(this, P)
            this.P = P;
        end
        
        % Get and Set parameters
        function SetParam(this, params, values)
            this.P = SetParam(this.P,params, values);
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
                error('no signal to plot. Use Sim command first.')
            end
      
            figure;
            h = SplotVar(this.P, varargin{:});
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
%           this.ResetParamSet();
            this.P = QuasiRefine(this.P,delta);
        end
    
        % Get the number of param vectors
        function nb_pts = GetNbParams(this)
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
 
    end
end