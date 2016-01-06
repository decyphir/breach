classdef BreachOpenSystem < BreachSystem
    % BreachOpenSystem is a BreachObject with an input generator 
    % TODO Documentation
       
    properties
        InputMap       % Maps input signals to idx in the input generator
        InputGenerator % BreachObject responsible for generating inputs
    end
    
    methods
        
        
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
                idx = this.InputMap(input{1});
                U.u(:,idx_mdl) = this.InputGenerator.P.traj.X(idx,:)';
            end
        end     
        
    end
    
end
