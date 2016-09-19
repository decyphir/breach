classdef DocRun
    % DocRun a small class to run, test, and publish scripts. 
    
    properties
        script_path
        script_dir
        script_name
        ref_result
   
        report
    end
    
    methods
        function this = DocRun(script, run_now)
            % DocRun(scriptname, run_now) -- if run_now is 1, runs as ref, otherwise run_and_cmp
            
            this.script_path = which(script);
            [this.script_dir, this.script_name] = fileparts(this.script_path);
            this.ref_result = [pwd filesep this.script_name '_ref_result.mat'];
           
            if nargin==2
                if run_now==1
                    this.run_as_ref();
                else
                    this= this.run_and_cmp();
                end
                
            end
            
            
        end
        
        function run(this)
            % run run the script as normal (that really useful?)            
            crd = pwd;
            cd(this.script_dir)
            evalin('base',this.script_name);
            cd(crd);
        end
        
        function run_as_ref(this)
            % run_as_ref run and save the resulting workspace
            crd = pwd;
            cd(this.script_dir)
            evalin('base',this.script_name);
        %    if input(['Saving ' this.ref_result ', are you sure (0 or 1)?\n'])
                fprintf('Saving...\n');
                evalin('base', ['save(''' this.ref_result ''');']);
        %    end
            cd(crd);
        end
        
        function this = run_and_cmp(this)
            %run_and_cmp run and compare to reference data
            
            % run script 
            crd = pwd;
            cd(this.script_dir);
            evalin('base',this.script_name);
   
            % Tries to load previous results
            try
                assignin('base','ref_result',load(this.ref_result));
            catch
                error('Reference results failed to load. Run run_as_ref to create.' )
            end

            % Compare with reference results
            fn = evalin('base','fieldnames(ref_result)');
            for i_f = 1:numel(fn)
                if evalin('base',['isa(ref_result.(''' fn{i_f} '''), ''BreachSet'')']);
                    cmp___ = evalin('base',[fn{i_f} '.compare(ref_result.(''' fn{i_f} '''))']);
                    this.report.(fn{i_f}) = cmp___;
                end
            end
            cd(crd)
        end
        
        function show_report(this)
           
            flds = fieldnames(this.report);
            for ifn = 1:numel(flds) 
                fn = this.report.(flds{ifn});
                sfn = sprintf([flds{ifn} ':']);
                ufn = regexprep(sfn,'.','-');
                fprintf([sfn '\n', ufn, '\n']);
                fn.printStatus;
                fprintf('\n');
            end
            
        end
    end
    
end