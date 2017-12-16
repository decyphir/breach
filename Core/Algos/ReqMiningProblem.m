classdef ReqMiningProblem < BreachCEGIS 
    
   methods
       function this = ReqMiningProblem(varargin)
           switch nargin
               case 2
                   % syntax: ReqMiningProblem(synth_pb, falsif_pb) or
                   % ReqMiningProblem(falsif_pb, synth_pb)
                   
                   pb1 = varargin{1}; 
                   pb2 = varargin{2};
                   
                   if isa(pb1, 'FalsificationProblem')&&isa(pb2, 'ParamSynthProblem')
                       this.synth_pb  = pb2;
                       this.falsif_pb = pb1;
                   elseif  isa(pb2, 'FalsificationProblem')&&isa(pb1, 'ParamSynthProblem')
                       this.synth_pb  = pb1;
                       this.falsif_pb = pb2;
                   elseif isa(pb1, 'BreachSystem')&&isa(pb2, 'STL_Formula')
                       Br = pb1;
                       phi = pb2;
                       this.falsif_pb = FalsificationProblem(Br, phi);
                       this.synth_pb  = ParamSynthProblem(Br, phi);
                   else
                       this.synth_pb  = pb1;
                       this.falsif_pb = pb2;
                       error('ReqMiningProblem:unknown_arg_type', 'Arguments should be of class FalsificationProblem and ParamSynthProblem.');
                   end
                   
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
        
       function new = copy(this)
            % copy operator, works with R2010b or newer.
            objByteArray = getByteStreamFromArray(this);
            new = getArrayFromByteStream(objByteArray);
        end

       
   end
end