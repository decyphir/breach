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
        report
    end
    
    methods
        %% Constructor
        function this = DocRun(script, run_now)
            % DocRun(scriptname, run_now) -- if run_now is 1, runs as ref, otherwise run_and_cmp
            InitBreach;
            this.current_dir = pwd;
            this.publish_dir = this.current_dir;
            
            this.script_path = which(script);
            % checks whether the script is in a namespace - e.g. BrDemo
            if (regexp(this.script_path,[filesep '\+']))
                % will run the script in the current folder
                this.script_dir = pwd;
                this.script_name = script;
            else
                % otherwise, script is run in its own folder
                [this.script_dir, this.script_name] = fileparts(this.script_path);
            end
            this.ref_result = [this.current_dir filesep this.script_name '_ref_result.mat'];
            
            if nargin==2
                if run_now==1
                    this.run_as_ref();
                else
                    this= this.run_and_cmp();
                end
            end
            
            this.report.status = BreachStatus();
            this.publish_src_dir = [this.current_dir filesep 'publish_src'];
            this.publish_html_dir = [this.current_dir filesep 'html'];
            
        end
        
        %% Run
        function result = run(this)
            % run run the script as normal, returns the output workspace result in a struct
            this.pre_run()
            evalin('base',this.script_name);
            evalin('base', 'save(''_tmp_docrun_save_'')');
            result = load('_tmp_docrun_save_');
            delete('_tmp_docrun_save_.mat');
            this.post_run();
        end
        
        function run_as_ref(this)
            % run_as_ref run and save the resulting workspace
            this.pre_run()
            evalin('base',this.script_name);
            fprintf('Saving %s...\n', this.ref_result);
            evalin('base', ['save(''' this.ref_result ''');']);
            this.post_run()
        end
             
        function run_and_cmp(this)
            
            % run script
            result= this.run();
            
            % Compare with reference results
            this.cmp_to_ref(result);
            
        end
        
        function cmp_to_ref(this, test_res)
            
            % Tries to load previous results
            try
                ref = load(this.ref_result);
            catch
                this.post_run();
                error('Reference results failed to load. Run run_as_ref to create.' )
            end
            
            fn_ref = fieldnames(ref);
            
            % go through results in ref, and make sure they are consistent
            for i_f = 1:numel(fn_ref)
                f = fn_ref{i_f};
                v_ref = ref.(f);
                if ~isfield(test_res,f)
                    this.report.status.addStatus(-1, ['Variable ' f ' missing']);
                elseif ~isequal(ref.(f), test_res.(f))
                    v_test = test_res.(f);
                    if any(strcmp('compare',methods(v_ref)))
                        cmp = v_ref.compare(v_test);
                        if (cmp.status ~= 0)
                            this.report.status.addStatus(-1,['Variable ' f ' is different']);
                            this.report.diffs.(f) = cmp;
                        end
                    else
                        this.report.status.addStatus(-1,['Variable ' f ' is different']);
                    end
                end
            end
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
            % saves workspace
            evalin('base','save docrun_backup.mat');
            evalin('base','clear');          
            cd(this.script_dir);
        end
        
        function [success, msg, msg_id] = post_run(this)
            % checks error status 
            success = 1; 
            [msg, msg_id] = lasterr;           
            
            % get back to current folder, and load back workspace
            cd(this.current_dir);
            evalin('base','load docrun_backup.mat');
            delete('docrun_backup.mat')
        end
        
        %% Publish
        function success = publish_beamer(this, compile, open)
            % publish_beamer run and create beamer
            success = 1;
            global BreachGlobOpt;
            this.pre_run()
            if ~(exist(this.publish_src_dir,'dir'))
                mkdir(this.publish_src_dir);
            end
            
            % get beamer template
            breach_publish_stuff_dir = [BreachGlobOpt.breach_dir filesep 'Core' filesep 'm_src' filesep 'publish_stuff' filesep];
               
            % run publish command
            res = evalin('base', ['publish(''' this.script_name ''',' 'struct(''format'',''latex'',''catchError'',false,''stylesheet'',''' breach_publish_stuff_dir  'matlab2beamer.xsl'', ''outputDir'',''' this.publish_src_dir ''')' ');']);
            
            % copy Decyphir Logo
            copyfile([breach_publish_stuff_dir filesep 'DecyphirLogo.png'], [this.publish_dir filesep 'publish_src']); 
            
            % get resulting file
            [tex_src_path, tex_src]= fileparts(res);
            pdf_file = [tex_src '.pdf'];
            tex_file = [tex_src '.tex'];
            
            % clean up
            this.post_run();
            
            % FIXME the following only works on m
            if nargin>1
                
                if ~(exist(this.publish_dir,'dir'))
                    mkdir(this.publish_dir);
                end
                
                if compile
                    cd(tex_src_path);
                    [status, result] = system(['pdflatex ' tex_file]);
                    if (status)
                        success = 0; 
                        warning('Couldn''t compile Beamer presentation.');
                    else
                        copyfile(pdf_file, [this.publish_dir filesep this.script_name '.pdf']);
                    end
                    cd(this.current_dir);
                end
                
                if nargin>2
                    if open
                        cd(this.publish_dir);
                        if ispc
                            system(['"c:\Program Files\SumatraPDF\SumatraPDF.exe" ' this.script_name '.pdf']);
                        else
                            system(['open ' this.script_name '.pdf']);
                        end
                        
                        cd(this.current_dir);
                    end
                end
                
            end
            
        end
        
        function publish_html(this, op)
            
            this.pre_run()
            if ~(exist(this.publish_html_dir,'dir'))
                mkdir(this.publish_html_dir);
            end
            
            % run publish command
            res = evalin('base', ['publish(''' this.script_name ''',' 'struct(''format'',''html'', ''catchError'',false,''outputDir'',''' this.publish_html_dir ''')' ');']);
            
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