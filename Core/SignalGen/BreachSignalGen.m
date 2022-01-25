classdef BreachSignalGen < BreachSystem
 % BreachSignalGen A class to generate signals of different types. 
 %   This class is derivated from BreachSystem, and thus inherits from all properties and methods.  
 %   It aggregates several instances of a simpler signal_gen class. The
 %   main use-case of this class is to as input generators for a
 %   BreachOpenSystem. It can also be used to interface an external
 %   simulator. 
 % 
 % BreachSignalGen Properties
 %   signalGenerators - cell array of signalgen object.
 %   dt_default=1e-3  - default fixed time step for signal generation.
 % 
 % BreachSignalGen Methods
 %       BreachSignalGen - constructor, takes a cell array of signal_gen
 %                         objects as only argument
 %
 %  
 % See also BreachOpenSystem, signal_gen
 
    properties
        signalGenerators
        dt_default=1e-3 % in case no time step is provided
    end
    
    methods
        %% Constructor
        function this = BreachSignalGen(signalGenerators)
    
            if nargin==0
               return; 
            end
            
            if ~iscell(signalGenerators)
               signalGenerators = {signalGenerators}; 
            end
            
            this.Domains = [];
            this.InitSignalGen(signalGenerators);
            this.CheckinDomain();
        end
        
        function InitSignalGen(this, signalGenerators)
            this.signalGenerators = signalGenerators;
            signals ={}; 
            params = {};
            p0=[];

            ParamDomains = []; 
            SignalDomains = []; 
            for isg = 1:numel(signalGenerators)
                sg=  signalGenerators{isg};
                signals = {signals{:}, signalGenerators{isg}.signals{:}};
                params = {params{:}, signalGenerators{isg}.params{:}}; 
                
                % domains 
                num_sig = numel(sg.signals);
                if isempty(sg.signals_domain)
                    SignalDomains = [SignalDomains repmat(BreachDomain(),1, num_sig)];
                else
                    SignalDomains = [SignalDomains sg.signals_domain];
                end
                
                num_par = numel(sg.params);
                if isempty(sg.params_domain)
                    ParamDomains = [ParamDomains repmat(BreachDomain(),1, num_par)];
                else
                    ParamDomains = [ParamDomains sg.params_domain];
                end
                
                % default values
                p0sg = signalGenerators{isg}.p0;
                if size(p0sg,2) >1
                    p0sg = p0sg';
                end
                p0 = [p0; p0sg ];
            end
            
            %% Uniquify
            [signals, ius] = unique(signals,'stable');
            SignalDomains =SignalDomains(ius);
            [params, iup] = unique(params,'stable');
            ParamDomains =ParamDomains(iup);
            
            
            this.Domains = [SignalDomains ParamDomains];
            
            p0 = [zeros(numel(signals),1) ; p0 ];
            this.Sys = CreateExternSystem('BreachSignalGen', signals, params, p0, @(Sys, tspan, p)breachSimWrapper(this, Sys, tspan, p));
            this.Sys.tspan =0:.01:10;
            this.Sys.Verbose =0;
            this.P = CreateParamSet(this.Sys);
            this.P.epsi(:,:)=0;
            
            if isaSys(this.Sys) % Note: we ignore initial conditions for now in ParamRanges
                                       % OK for Simulink, less so for ODEs...
                this.SignalRanges = [];
            end
            
            
        end
        
        function [tspan, X] = breachSimWrapper(this, Sys, tspan, p)
            
            if numel(tspan)==1
               tspan = 0:this.dt_default:tspan; 
            elseif numel(tspan)==2
               tspan = tspan(1):this.dt_default:tspan(2); 
            end
            
            cur_is =1;
            
            for isg = 1:numel(this.signalGenerators)
               sg = this.signalGenerators{isg};
               idx_psg = FindParam(this.P, sg.params);
               p_isg = p(idx_psg);
               ns = numel(sg.signals);
               [X(cur_is:cur_is+ns-1, :), tspan]= sg.computeSignals(p_isg, tspan);
               cur_is = cur_is+ ns;
            end
                        
        end
        
        function SetParam(this, params, values, varargin)
            SetParam@BreachSet(this, params, values, varargin{:});
            pts = [];
            
            % look for parameters from_file_signal generators
            for isg = 1:numel(this.signalGenerators)
                sg = this.signalGenerators{isg};
                if isa(sg, 'from_file_signal_gen')
                    % parameters in from_file_signal_gen are enum, listed
                    % in a field pts. First element of params is going to
                    % determine the others
                    if iscell(params)
                        param_ref = params{1};
                    else
                        param_ref = params;
                    end
                    if isnumeric(param_ref)
                        param_ref = this.P.ParamList{param_ref};
                    end
                    pts = GetParam(this, sg.params);
                    
                    i_ref = find(strcmp(param_ref,sg.params),1);
                    if ~isempty(i_ref) % might be a requirement parameter
                        for jp = 1:size(pts,2)
                            p = pts(:, jp);
                            p_is_good= false;  % checks whether p is a valid parameter, i.e., can be found in sg.pts
                            for sg_jp = 1:size(sg.pts,2) % look into columns of signal_generator matrix of valid parameters
                                if isequal(p, sg.pts(:,sg_jp));
                                    p_is_good = true;
                                    break;
                                end
                            end
                            if p_is_good
                                continue;
                            else  % we need to fix p. First find a valid parameter value for param_ref
                                pref =p(i_ref);
                                j_ref = find(sg.pts(i_ref,:)==pref,1);
                                if isempty(j_ref) % find closest - should be the job of checkin domain but whatever... 
                                    [~, j_ref] = min(abs(sg.pts(i_ref,:)-pref));
                                    p(i_ref) = sg.pts(i_ref, j_ref); 
                                end
                                p = sg.pts(:, j_ref);
                                pts(:, jp) = p;
                            end
                        end
                        this.SetParam@BreachSet(sg.params, pts);
                    end
                end
            end
            
        end
        
        function sg = GetSignalGenFromSignalName(this, sig_name)
            % GetSignalGenFromSignalName(this, sig_name) only works for one
            % dimensional signal generators so far
            for  is = 1:numel(this.Sys.DimX)
                for isg = 1:numel(this.signalGenerators)
                    ll = strcmp(this.signalGenerators{isg}.signals, sig_name);
                    if any(ll)
                        sg = this.signalGenerators{isg};
                        % fix domains
                        sg.params_domain = this.GetDomain(sg.params);
                        sg.signals_domain = this.GetDomain(sg.signals);                        
                        return;
                    end
                end
                error('GetSignalGenFromSignalName:signal_not_found.', 'No signal generator generates this signal.');
            end
        end
        
        function signal_gen_cfg = ExportConfig(this)

            signal_gen_cfg = {};
            
            for isg = 1:numel(this.signalGenerators)
                sg =this.signalGenerators{isg};
                cfg = struct;
                cfg.signal_gen_class = class(sg);
                args_name = sg.getSignalGenArgs();
                cfg.Args = {sg.signals};
                for ia= 1:numel(args_name)
                    cfg.Args{end+1} = sg.(args_name{ia});
                end
                if size(sg.p0,2)>1
                    cfg.Args{end+1} = sg.p0;
                else
                    cfg.Args{end+1} = sg.p0';
                end
                signal_gen_cfg{end+1} = cfg; 
            end
            
        end
        
        
    end
end
