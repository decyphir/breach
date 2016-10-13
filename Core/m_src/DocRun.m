classdef DocRun
    % DocRun a class to run, test, and publish scripts. Still quite
    % hack-ish. Allows to run things in relatively clean/controlled
    % workspace to improve chance or repeatability (limit side-effects).
    % Can save results and compare, and publish beamer slides.
    
    properties
        script_path
        script_dir
        script_name
        ref_result
        current_dir  %  dir where the object is created
        report
    end
    
    methods
        function this = DocRun(script, run_now)
            % DocRun(scriptname, run_now) -- if run_now is 1, runs as ref, otherwise run_and_cmp
            
            this.current_dir = pwd;
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
            this.pre_run()
            evalin('base',this.script_name);
            this.post_run();
        end
        
        function run_as_ref(this)
            % run_as_ref run and save the resulting workspace
            crd = pwd;
            cd(this.script_dir)
            evalin('base','save docrun_backup.mat');
            
            evalin('base',this.script_name);
            if input(['Saving ' this.ref_result ', are you sure (0 or 1)?\n'])
                fprintf('Saving...\n');
                evalin('base', ['save(''' this.ref_result ''');']);
            end
            cd(crd);
        end
        
        function this = run_and_cmp(this)
            %run_and_cmp run and compare to reference data
            
            % run script
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
            evalin('base','load docrun_backup.mat');
            delete('docrun_backup.mat')
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
        
        function publish_beamer(this, compile, open)
            global BreachGlobOpt;
            this.pre_run()
            publish_tmp_dir = [this.current_dir filesep 'publish_tmp' filesep];
            if ~(exist(publish_tmp_dir,'dir'))
                mkdir(publish_tmp_dir);
            end
            breach_publish_stuff_dir = [BreachGlobOpt.breach_dir filesep 'Core' filesep 'm_src' filesep 'publish_stuff' filesep];
            evalin('base', ['publish(''' this.script_name ''',' 'struct(''format'',''latex'',''stylesheet'',''' breach_publish_stuff_dir  'matlab2beamer.xsl'', ''outputDir'',''' publish_tmp_dir ''')' ')']);
            this.post_run()
            if nargin>1
                if compile
                    cd(publish_tmp_dir);
                    system(['pdflatex ' this.script_name '.tex']);
                    cd('..');
                end
                if nargin>2
                    if open
                        cd(publish_tmp_dir);
                        system(['open ' publish_tmp_dir filesep this.script_name '.pdf']);
                        cd('..');
                    end
                end
            end
            
            
        end
        
        function compile_beamer(this)
            publish_tmp_dir = [this.current_dir filesep 'publish_tmp' filesep];
            cd(publish_tmp_dir);
            system(['pdflatex ' this.script_name '.tex']);
            cd('..');
        end
        
        function open_beamer(this)
            publish_tmp_dir = [this.current_dir filesep 'publish_tmp' filesep];
            system(['open ' publish_tmp_dir filesep this.script_name '.pdf']);
            cd('..');
        end
        
        function pre_run(this)
            evalin('base','save docrun_backup.mat');
            cd(this.script_dir);
        end
        
        function post_run(this)
            cd(this.current_dir);
            evalin('base','load docrun_backup.mat');
            delete('docrun_backup.mat')
        end
        
    end
    
end