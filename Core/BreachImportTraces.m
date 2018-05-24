classdef BreachImportTraces < BreachTraceSystem
    
    methods
        function this =BreachImportTraces(traces)
            
              if ~exist('traces', 'var')||isempty(traces)
                [filenames, paths] = uigetfile( ...
                    {  '*.mat','MAT-files (*.mat)'}, ...
                    'Pick one or more files', ...
                    'MultiSelect', 'on');
                
                if isequal(filenames,0) % cancel
                    fname = {};
                else
                    if ~iscell(filenames)
                        filenames= {filenames};
                    end
                    fname = cellfun( @(c)([ paths c  ] ), filenames,'UniformOutput',false);
                end
              end
              for it = 1:numel(traces)
                    traces{it} = matfile(traces{it}, 'Writable', true);
              end
       
              if isstruct(traces)
                  traces =  {traces};
              end
              
              signature = traces{1}.signature;
              this.Sys = CreateExternSystem('TraceObject', signature.signals,  signature.params);
              this.P = CreateParamSet(this.Sys);
 
                %  Default domains
                for i = 1:this.P.DimX
                this.Domains(i) = BreachDomain(signature.signal_types(i));
                end
                for i = this.P.DimX+1:this.P.DimP
                    this.Domains(i) = BreachDomain(signature.param_types(i-this.P.DimX));
                end
          
                for it = 1:numel(traces)
                    this.AddTrace(traces{it});
                end
        end
    
        function AddTrace(this, traj)
                    
            Pnew = CreateParamSet(this.Sys);
            Pnew.epsi(:,:) = 0;
            
            nb_traces =this.CountTraces();
            
            Pnew.traj={traj};
            Pnew.traj_ref = 1;
            Pnew.traj_to_compute =  [];
            Pnew.pts(1:Pnew.DimP,1) =  traj.param';
            if nb_traces == 0
                this.P = Pnew;
            else
                this.P = SConcat(this.P, Pnew);
            end
            this.P.traj_ref = 1:nb_traces+1;
            this.P.traj_to_compute =  [];
            this.Sys.tspan = traj.time;
            
            if isfield(this.P, 'Xf')
                this.P.Xf(:,end+1)= traj.X(:,end);
            else
                this.P.Xf= traj.X(:,end);
            end
            
        end
    end
    
    
end