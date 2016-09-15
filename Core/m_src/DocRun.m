classdef DocRun
    
    properties
        script_path
        script_name
        data_ref_path
    end
    
    methods
        function this = DocRun(script)
            this.script_path = which(script);
            [~, this.script_name] = fileparts(this.script_path);
            this.data_ref_path = regexprep(this.script_path,'.m$', '_data_ref.mat');
        end
        
        
        function run(this)
            % run run the script as normal (that really useful?)
            run(this.script_path);
        end
        
        function run_as_ref(this)
            % run_as_ref run and save the resulting workspace
            run(this.script_path);
            save(this.data_ref_path);
        end
        
        function report = run_and_cmp(this)
            % run and compare to reference data
            
            try
                data_ref = load(this.data_ref_path);
            catch
                error('Feference data failed to load. Run run_as_ref to create.' )
            end
            
            run(this.script_path);
            
            fn = fieldnames(data_ref);
            for i_f = 1:numel(fn)
                
                if isa(data_ref.(fn{i_f}), 'BreachSystem')
                    eval(['cmp___=' fn{i_f} '.compare(data_ref.(''' fn{i_f} '''));']);
                    report.(fn{i_f}) = cmp___;
                end
            end
        end
        
    end
    
end