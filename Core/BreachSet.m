classdef BreachSet < handle
    properties
        P % legacy Breach structure, contains data and trajectories, etc
    end
    
    methods
        
        % Get and Set of the parameter set structure
        function P = GetP(this)
            P = this.P;
        end
        function SetP(this, P)
            this.P = P;
        end
        
               % Get and Set parameters for simulation
        function SetParam(this, params, values)
            this.P = SetParam(this.P,params, values);
            this.UpdateParamRanges();
        end
        function values = GetParam(this, params)
            values = GetParam(this.P,params);
        end
        
        % Set param ranges
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
        
        % Grid Sample - reset P (?)
        function GridSample(this, delta)
%           this.ResetParamSet();
            this.P = Refine(this.P,delta);
        end
        
        % Quasi-Random Sample - reset P
        function QuasiRandomSample(this, delta)
%           this.ResetParamSet();
            this.P = QuasiRefine(this.P,delta);
        end
    
        % Get the number of param vectors
        function nb_pts=  GetNbParams(this)
            nb_pts= size(this.P.pts,2);
        end
        
        
        % Plot parameters
        function PlotParams(this, varargin)
            figure;
            SplotPts(this.P, varargin{:});
        end
 
    end
end