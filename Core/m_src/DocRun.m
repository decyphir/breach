classdef DocRun < handle
    % DocRun a class to run, test, and publish scripts. Still quite
    % hack-ish. Allows to run things in relatively clean/controlled
    % workspace to improve chance or repeatability (limit side-effects).
    % Can save results and compare, and publish beamer slides.
    
    properties
        script_name
        script_path
        script_dir
        ref_result
        current_dir  %  dir where the object is created
        publish_dir
        publish_src_dir
        publish_html_dir
    end
    
    methods
        %% Constructor
        function this = DocRun(script, run_now)
            % DocRun(scriptname, run_now) -- if run_now is 1, runs as ref, otherwise run_and_cmp
            
            this.current_dir = pwd;
            this.publish_dir = this.current_dir;
            
            this.script_path = which(script);
            % checks whether the script is in a namespace
            if (regexp(this.script_path,[filesep '\+']))
                % will run the script in the current folder
                this.script_dir = pwd;
                this.script_name = script;
            else
                % otherwise, script is run in its own folder
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
            
            this.publish_src_dir = [this.current_dir filesep 'publish_src'];
            this.publish_html_dir = [this.current_dir filesep 'html'];
            
        end
        %% Run
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
            % shows comparison results
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
        
        function pre_run(this)
            % ensures
            evalin('base','save docrun_backup.mat');
            cd(this.script_dir);
        end
        
        function post_run(this)
            cd(this.current_dir);
            evalin('base','load docrun_backup.mat');
            delete('docrun_backup.mat')
        end
        
        %% Publish
        function publish_beamer(this, compile, open)
            % publish_beamer run and create beamer
            global BreachGlobOpt;
            this.pre_run()
            if ~(exist(this.publish_src_dir,'dir'))
                mkdir(this.publish_src_dir);
            end
            
            % get beamer template
            breach_publish_stuff_dir = [BreachGlobOpt.breach_dir filesep 'Core' filesep 'm_src' filesep 'publish_stuff' filesep];
            
            % run publish command
            res = evalin('base', ['publish(''' this.script_name ''',' 'struct(''format'',''latex'',''stylesheet'',''' breach_publish_stuff_dir  'matlab2beamer.xsl'', ''outputDir'',''' this.publish_src_dir ''')' ');']);
            
            % get resulting file
            [tex_src_path, tex_src]= fileparts(res);
            pdf_file = [tex_src '.pdf'];
            tex_file = [tex_src '.tex'];
            
            % clean up
            this.post_run()
            
            % FIXME the following only works on macOS
            if nargin>1
                
                if ~(exist(this.publish_dir,'dir'))
                    mkdir(this.publish_dir);
                end
                
                if compile
                    cd(tex_src_path);
                    [status] = system(['pdflatex ' tex_file]);
                    if (status)
                        error('Couldn''t compile Beamer presentation.');
                    else
                        copyfile(pdf_file, [this.publish_dir filesep this.script_name '.pdf']);
                    end
                    cd(this.current_dir);
                end
                if nargin>2
                    if open
                        cd(this.publish_dir);
                        system(['open ' this.script_name '.pdf']);
                        cd(this.current_dir);
                    end
                end
            end
            
        end
        
        function publish_html(this, op)
            global BreachGlobOpt;
            this.pre_run()
            if ~(exist(this.publish_html_dir,'dir'))
                mkdir(this.publish_html_dir);
            end
            
            % run publish command
            res = evalin('base', ['publish(''' this.script_name ''',' 'struct(''format'',''html'', ''outputDir'',''' this.publish_html_dir ''')' ');']);
            
            % open
            if exist('op','var')
                if op
                    open(res);
                end
            end
            
            % clean up
            this.post_run();
        end
        
        %% For release
        
        function set_release_dir(this)
            
            global BreachGlobOpt;
            this.publish_html_dir = [BreachGlobOpt.breach_dir filesep 'Doc' filesep 'html'];
            this.publish_dir = [BreachGlobOpt.breach_dir filesep 'Doc' filesep 'pdf'];
            
        end
    end
end