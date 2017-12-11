classdef BreachImportData < BreachSignalGen
% BreachImportData Specialisation of BreachSignalGen to
% from_file_signal_gen - main thing is that when importing signals, it can use
% different times, also pre-sample and loads all traces 
      
properties
    fname
end

methods
    function this = BreachImportData(fname, signals) 

        if ~exist('fname', 'var')||isempty(fname)
            
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
        
        if ~exist('signals', 'var')
            signals = {};
        end  
        
        if ~isempty(fname)
            ff_ingen = from_file_signal_gen(signals, fname);
        else
            ff_ingen = constant_signal_gen('x');
        end
        
        this = this@BreachSignalGen({ff_ingen});
        
        if ~isa(ff_ingen, 'constant_signal_gen')  % not canceled 
            dom = this.GetDomain('file_idx');  % will need to improve this at some point
            dom.domain = [dom.enum(1) dom.enum(end)];
            this.SetDomain('file_idx', dom);
            this.SampleDomain('file_idx', 'all');
            this.Sys.Verbose=0;
            this.Sim();
        end

        this.fname = fname;
    end
    
    function [tspan, X] = breachSimWrapper(this, Sys, tspan, p)
            
             sg = this.signalGenerators{1};
 
            if numel(tspan)==1
               tspan = 0:this.dt_default:tspan; 
            elseif numel(tspan)==2
               tspan = tspan(1):this.dt_default:tspan(2); 
            end
            
            % Needs some more cleanup
            p = p(this.Sys.DimX+1:end);
            cur_ip =1;
            cur_is =1;
            
            np = numel(sg.params);
            p_isg = p(cur_ip:cur_ip+np-1);  %
            ns = numel(sg.signals);
            [X(cur_is:cur_is+ns-1, :), tspan] = sg.computeSignals(p_isg, tspan);
            
    end
    
            function [summary, traces] = ExportTracesToStruct(this,i_traces, varargin)
            % BreachImportData.ExportTracesToStruct
            
            summary = [];
            traces = [];
            if ~this.hasTraj()
                error('Breach:ExportTrace:no_trace', 'No trace to export - run Sim command first');
                return;
            end
            
            num_traces = numel(this.P.traj);
            if nargin==1
                i_traces = 1:num_traces;
            end
            
            % Additional options
            options = struct('FolderName', []);
            options = varargin2struct(options, varargin{:});
            
            if isempty(options.FolderName)
                options.FolderName = ['Import_Results_' datestr(now, 'dd_mm_yyyy_HHMM')];
            end
            
            %% Common stuff
            
            % parameter names
            param_names = this.GetSysParamList();
            
            % input signal names
            signal_names = this.GetSignalList();
                      
            if isfield(this.P,'props_names')
                spec_names = this.P.props_names;
            end

            summary.date = datestr(now);
            summary.num_traces = num_traces;
            summary.params.names = param_names;
            summary.params.values = this.GetParam(summary.params.names);

            
            %% traces
            file_idx= this.GetParam('file_idx');
            summary.filenames = {};
            summary.paths = {};
            for it = i_traces
                
                % filename 
                ip = file_idx(it);
                [pth, fn] = fileparts(this.signalGenerators{1}.file_list(ip).name); 
                traces(it).filename = fn;
                traces(it).path = pth;
                
                summary.filenames   = [summary.filenames traces(it).filename];
                summary.paths   = [summary.paths traces(it).path];
                % params
                traces(it).params.names = param_names;
                traces(it).params.values = this.GetParam(param_names,it)';
                
                % time
                traces(it).time = this.P.traj{it}.time;
                
                % input signals
                traces(it).signals.names = signal_names;
                traces(it).signals.values =  this.GetSignalValues(signal_names, it);
                                
                % specifications
                if isfield(this.P,'props_names')
                    traces(it).specs.ids = spec_names;
                    for ip = 1:numel(this.P.props_names)
                        traces(it).specs.stl_formula{ip} = disp(this.P.props(ip));
                        traces(it).specs.stl_formula_full{ip} = disp(this.P.props(ip),0);
                        params = get_params(this.P.props(ip));
                        traces(it).specs.params(ip).names = fieldnames(params);
                        traces(it).specs.params(ip).values = this.GetParam(fieldnames(params), it)';
                        traces(it).specs.rob(ip).time =this.P.props_values(ip, it).tau;
                        traces(it).specs.rob(ip).values =  this.P.props_values(ip, it).val;
                        traces(it).specs.status(ip) =  this.P.props_values(ip, it).val(1)>=0;
                    end
                end
            end
            
            if isfield(this.P, 'props')
                summary.specs.names = spec_names;
                this.SortbyRob();
                this.SortbySat();
                summary.specs.rob = this.GetSatValues();
                summary.specs.sat = summary.specs.rob>=0;
                summary.num_sat = - sum( ~summary.specs.sat, 1  );
            end
        end
        
        function ExportToExcel(this, varargin)
            
            % Additional options
            options = struct('FileName', 'Results.xlsx');
            options = varargin2struct(options, varargin{:});
            
            global BreachGlobOpt
            
            summary = this.ExportTracesToStruct();
            
            if ~isfield(summary, 'specs')
                warning('Breach:ExportToExcel:no_spec_for_Excel','Export to Excel requested but there is no requirement result to report. Excel file not created.');
            else
                excel_file =options.FileName;
                breach_dir = BreachGlobOpt.breach_dir;
                template_file_path = [breach_dir filesep 'Ext' filesep 'Toolboxes' filesep 'ExportResults' filesep 'BreachResults_template.xlsx'];
                copyfile(template_file_path, excel_file);
                
                % Write header
               
                for ispec = 1:numel(summary.specs.names)
                    hdr{ispec} = ['Req. ' num2str(ispec)];
                end
                
                for iparam = ispec+1:ispec+numel(summary.params.names)
                    hdr{iparam} = ['param. ' num2str(iparam-ispec) ];
                end
                
                xlswrite(excel_file, hdr, 1, 'D1');
                xlswrite(excel_file, ['Path' 'Filename' summary.specs.names summary.params.names], 1, 'B2');
                 
                xlswrite(excel_file, summary.paths', 1, 'B3');
                xlswrite(excel_file, summary.filenames', 1, 'C3');
                xlswrite(excel_file,  summary.num_sat'  , 1, 'A3');
  
                % Write data
                xlswrite(excel_file, [ summary.specs.rob' summary.params.values'] , 1, 'D3');
                
                this.disp_msg(['Summary written into ' excel_file]);
            end
            
        end
        
    
end
end