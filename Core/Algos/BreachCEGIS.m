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
            while (cont)
                clc;
                %% Synthesis step
                fprintf('Iter %d/%d\n', this.iter,this.iter_max)
                fprintf('Synthesis step\n');
                fprintf('--------------\n');
                this.synth_pb.solve();
                
                if isempty(this.synth_pb.obj_best<0)
                    fprintf('Couldn''t synthesize parameters.\n');
                    return
                end
                
                %% Falsification step
                fprintf('Counter-Example step\n');
                fprintf('--------------------\n');
                this.falsif_pb.BrSet.SetParam(this.synth_pb.params, this.synth_pb.x_best, true);
                this.falsif_pb.ResetObjective();
                this.falsif_pb.solve();
                BrFalse = this.falsif_pb.GetBrSet_False();
                if isempty(BrFalse)
                    return
                end
                
                %% Update parameter synthesis problem
                this.synth_pb.BrSet.Concat(BrFalse.BrSet);
                this.synth_pb.ResetObjective();
                this.synth_pb.BrSys.Sys.Verbose=0;
                
                this.iter = this.iter+1;
                cont = this.iter<this.iter_max;
            end
        end
        
        function ResetIter(this)
           this.iter=0; 
        end
        
    end
    
end