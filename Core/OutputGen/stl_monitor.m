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
                this.formula= STL_Formula(formula);
            elseif isa(formula,'STL_Formula')
                this.formula= formula;    
            else
                error('stl_monitor:bad_constructor_arg', 'stl_monitor constructor requires a string or STL_Formula as argument.')                  
            end                    
            
            % collect signals and params names
            [this.in_signals, this.in_params, this.p0] = STL_ExtractSignals(this.formula);
          
            % construct legacy structures
            this.Sys = CreateExternSystem([this.formula_id '_Sys'], this.signals, {}, []);
            this.P = CreateParamSet(this.Sys);
           
            traj.param = zeros(numel(this.in_signals));
            traj.time = [];
            traj.X = [];
            traj.status = 0;
            
            this.P.traj = {traj};
            this.P.traj_ref = 1;
            this.P.traj_to_compute = [];
           
           
            % outputs
            this.out_signals = {get_id(this.formula)};
            this.out_values =  {[get_id(this.formula) '_rob']};        
            
            
        end
        
        function [rob, tau, val] = eval(this, time, X, p)
            this.P.traj{1}.X = X;
            this.P.traj{1}.time = time;
            [val, tau] = STL_Eval(this.Sys,this.P,this.P.traj{1}, time);  
        end
   end
end