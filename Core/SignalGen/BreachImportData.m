classdef BreachImportData < BreachSignalGen
% BreachImportData Specialisation of BreachSignalGen to
% from_file_signal_gen - main thing is that when importing signals, it can use
% different times, also pre-sample and loads all traces 
      
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
end
end