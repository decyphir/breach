classdef BreachCEGIS < BreachStatus
    
    properties
        synth_pb
        falsif_pb
        iter = 0
        iter_max = 10
    end
    
    methods
        function this = BreachCEGIS(SynthPb, FalsifPb)
            %BreachCEGIS class implements a Counter-Example Guided Inductive Synthesis strategy
            if nargin==0
                return;
            end
            this.synth_pb = SynthPb;
            this.falsif_pb = FalsifPb;
        end
        
        
        function solve(this)
            cont = true;
            BrFalse = []; 
            %% if Falsification not already done, start with that 
            if this.falsif_pb.nb_obj_eval ==0
                this.falsif_pb.solve();
                BrFalse = this.falsif_pb.GetFalse();
                if ~isempty(BrFalse)                    
                    this.synth_pb.ResetObjective(BrFalse.BrSet,this.synth_pb.params, [this.synth_pb.lb this.synth_pb.ub]);
                else
                    cont = false;
                end
            end
            %% Main Loop
            while (cont)
                clc;
                fprintf('Iter %d/%d\n', this.iter,this.iter_max)
                %% Synthesis step
                synth_step();
                
                %% Falsification step
                counter_ex_step();
                 
                %% Update parameter synthesis problem
                update_synth();
                this.iter = this.iter+1;
                cont = cont&&(this.iter<this.iter_max);
            end
            
            function update_synth()
                if cont
                    this.synth_pb.BrSet.Concat(BrFalse.BrSet);
                    this.synth_pb.ResetObjective();
                    this.synth_pb.BrSys.Sys.Verbose=0;
                end
            end
            
            function counter_ex_step()
                fprintf('Counter-Example step\n');
                fprintf('--------------------\n');
                this.falsif_pb.BrSet.SetParam(this.synth_pb.params, this.synth_pb.x_best, true);
                this.falsif_pb.ResetObjective();
                this.falsif_pb.solve();
                BrFalse = this.falsif_pb.GetFalse();
                if isempty(BrFalse)
                    cont = false;
                    return
                end
            end
            
            function synth_step()
                fprintf('Synthesis step\n');
                fprintf('--------------\n');
                this.synth_pb.solve();
                
                if isempty(this.synth_pb.obj_best<0)
                    fprintf('Couldn''t synthesize parameters.\n');
                    return
                end
                
            end
            
            
        end
        
        function ResetIter(this)
           this.iter=0; 
        end
        
    end
    
end