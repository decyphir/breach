classdef stl_monitor < output_gen

    properties
        P
        Sys
        formula
        formula_id
    end
    methods
        function this = stl_monitor(formula)
            if ischar(formula)
                this.formula= STL_Formula(STL_NewID('phi'), formula);
            elseif isa(formula,'STL_Formula')
                this.formula= formula;    
            else
                error('stl_monitor:bad_constructor_arg', 'stl_monitor constructor requires a string or STL_Formula as argument.')                  
            end                    
            
            % collect signals and params names
            [this.in_signals, this.in_params, this.p0] = STL_ExtractSignals(this.formula);
          
            % construct legacy structures
            this.Sys = CreateExternSystem([this.formula_id '_Sys'], this.in_signals, this.in_params, this.p0);
            this.P = CreateParamSet(this.Sys);
           
            traj.param = zeros(1,numel(this.in_signals)+numel(this.in_params));
            traj.time = [];
            traj.X = [];
            traj.status = 0;
            
            this.P.traj = {traj};
            this.P.traj_ref = 1;
            this.P.traj_to_compute = [];
           
            % outputs
            this.out_signals = {get_id(this.formula)};
            this.out_values =  {[get_id(this.formula) '_rob']};        
           
            % Init domains
            for vv =  [this.in_signals this.out_signals this.in_params this.out_values ]
                this.domains(vv{1}) = BreachDomain();
            end
            
        end
        
        function [rob, tau, val] = eval(this, time, X, p)
            this.P.traj{1}.X = X;
            this.P.traj{1}.time = time;
            if nargin>=4&&~isempty(p)
                this.P = SetParam(this.P, this.in_params,p);
            end
            [val, tau] = STL_Eval(this.Sys, this.formula, this.P,this.P.traj{1}, time);  
            rob = val(1);
        end
   end
end