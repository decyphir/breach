classdef BreachOpenSystem < BreachSystem
    % BreachOpenSystem is a BreachObject with an input generator
    % TODO Documentation
    
    properties
        InputMap       % Maps input signals to idx in the input generator
        InputGenerator % BreachObject responsible for generating inputs
    end
    
    methods
        
        function Sim(this,tspan,U)
            
            if ~exist('tspan','var')
                tspan = this.Sys.tspan;
            end
            Sys = this.Sys;
            if exist('U','var')
                if isnumeric(U)
                    DimU = this.InputMap.Count();
                    if size(U, 2)~=DimU+1;
                        err_msg= fprintf('Input must be an array with %d columns, first one being time.',DimU);
                        error(err_msg);
                    end
                    Us.t = U(:,1);
                    Us.u = U(:,2:end);
                else
                    Us = U;
                end
                Sys = this.Sys;
                Sys.init_u = @(~, pts, tspan) (Us);
            end
            
            this.P = ComputeTraj(Sys, this.P, tspan);
            
        end
           
        % we merge parameters of the input generator with those of the
        % system, but keep both BreachObjects
        function SetInputGen(this, IG)
            
            % First remove parameters from previous input generator
            Sys = this.Sys;
            idx_prev_inputs = this.GetParamsInputIdx();
            if (~isempty(idx_prev_inputs))
                idx_not_prev_inputs = boolean(ones(1, numel(Sys.ParamList)));
                idx_not_prev_inputs(idx_prev_inputs) = 0;
                Sys.ParamList = Sys.ParamList(idx_not_prev_inputs);
                Sys.p = Sys.p(idx_not_prev_inputs);
                Sys.DimP = numel(Sys.ParamList);
                this.Sys = Sys;
            end
            
            % IG can be a struct, a signal generator, or a BreachSystem
            
            %
            
            DimU = this.InputMap.Count;
            if (isstruct(IG))
                inputs = this.Sys.InputList;
                if ~isfield(IG,'type')
                    error('Input generator must be a struct with fields ''type'', ''cp'' at least')
                end
                
                if isscalar(IG.cp)
                    IG.cp = IG.cp*ones(1, DimU);
                end
                
                if ~isfield(IG,'method')
                    IG.method= 'previous';
                end
                
                if ischar(IG.method)
                    IG.method = {IG.method};
                end
                
                if numel(IG.method)==1
                    method = IG.method{1};
                    IG.method = cell(1,DimU);
                    for iu = 1:DimU
                        IG.method{iu} = method;
                    end
                end
                
                switch(IG.type)
                    case 'UniStep'
                        sg = fixed_cp_signal_gen(inputs, IG.cp, IG.method);
                        IG = BreachSignalGen({sg});
                    case 'VarStep'
                        sg = var_cp_signal_gen(inputs, IG.cp, IG.method);
                        IG = BreachSignalGen({sg});
                end
            else
                mm = methods(IG);
                if any(strcmp(mm,'computeSignals'))
                    IG = BreachSignalGen({IG});
                else
                    if ~any(strcmp(mm, 'Sim'))
                        error('Input generator should be a struct, a signal_gen derived object or a BreachSignalGen object');
                    end
                end
                
            end
            
            % Check Consistency - IG must construct signals for all signals in this.InputList
            for input = this.InputMap.keys
                idx= FindParam(IG.Sys,input);
                if idx<=IG.Sys.DimX
                    this.InputMap(input{1}) = idx;
                else
                    error(['Input ' input{1} ' is not provided by input generator.']);
                end
            end
            
            this.InputGenerator = IG;
            
            % Adds parameters for new input generator
            i_params = IG.Sys.DimX+1:IG.Sys.DimP;
            this.Sys = SetParam(this.Sys, IG.Sys.ParamList(i_params), IG.P.pts(i_params,1));
            this.Sys.DimP = numel(this.Sys.ParamList);
            
            % Resets P and ranges
            this.ParamRanges = [this.Sys.p this.Sys.p];
            this.SignalRanges = [];
            this.P = CreateParamSet(this.Sys);
            this.P.epsi(:,:) = 0;
            
            % Sets the new input function for ComputeTraj
            % FIXME?: tilde?
            this.Sys.init_u = @(~, pts, tspan) (InitU(this,pts,tspan));
            
        end
        
        function idx = GetParamsInputIdx(this)
            if isempty(this.InputGenerator)
                [~, idx] = FindParamsInput(this.Sys);
            else
                ig_params = this.InputGenerator.Sys.DimX+1:this.InputGenerator.Sys.DimP;
                idx = FindParam(this.Sys, this.InputGenerator.Sys.ParamList(ig_params));
            end
        end
        
        % calling the Input generator -
        % there might be saving to do if inputs are pre-generated
        function U = InitU(this, pts, tspan)
            
            idx_u = this.GetParamsInputIdx();
            ig_params = this.InputGenerator.Sys.DimX+1:this.InputGenerator.Sys.DimP;
            this.InputGenerator.P = SetParam(this.InputGenerator.P,ig_params,pts(idx_u));
            this.InputGenerator.Sim(tspan);
            U.t = this.InputGenerator.P.traj.time';
            U.u = zeros(numel(U.t), this.InputGenerator.Sys.DimX);
            idx_mdl = 0;
            for input= this.Sys.InputList % Sys.InputList is in the same order as the model
                idx_mdl = idx_mdl+1;
                idx =  FindParam(this.InputGenerator.P, input{1});
                U.u(:,idx_mdl) = this.InputGenerator.P.traj.X(idx,:)';
            end
        end
        
    end
    
end


