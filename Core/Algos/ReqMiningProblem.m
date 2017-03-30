classdef ReqMiningProblem < BreachCEGIS 
    
   methods
       function this = ReqMiningProblem(varargin)
           switch nargin
               case 2
                   % syntax: ReqMiningProblem(synth_pb, falsif_pb) 
                   this.synth_pb  = varargin{1};
                   this.falsif_pb = varargin{2};
                   
               case 4
                   % syntax: ReqMiningProblem(BrSet, phi,input_params, prop_params)
                   Br = varargin{1};
                   Br.Sim();
                   phi = varargin{2};
                   input_params= varargin{3};
                   prop_params= varargin{4};
                    
                   this.falsif_pb = FalsificationProblem(Br, phi, input_params.names, input_params.ranges);
                   this.synth_pb  = ParamSynthProblem(Br, phi, prop_params.names, prop_params.ranges);
                   this.synth_pb.setup_solver('binsearch');
  
           end
           rfprintf_reset();
       end
   
       function phi = GetRequirement(this)
          phi = this.synth_pb.Spec;
          p_mined =  this.synth_pb.x_best;
          phi = set_params(phi, this.synth_pb.params, p_mined);
       end
        
   end
end