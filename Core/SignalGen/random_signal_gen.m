classdef random_signal_gen < signal_gen
    % random_signal_gen  A class derived from signal_gen to generate signals from random control points
    %                    with variable time steps, using different interpolation methods.
    %
    % random_signal_gen Properties
    %   method - interpolation method for each signal, accepts all methods provided by interp1
    %
    % random_signal_gen Methods
    %   random_signal_gen -  constructor, takes signal names, number of
    %                        interpolation methods and a p0
    %  See also signal_gen.
    
    
    properties
        method   % interpolation method for each signal
    end
    
    methods
        
        function this = random_signal_gen(signals, method, p0)
            
            if nargin == 1
                method = 'previous';
            end
            
            if ~iscell(signals)
                signals = {signals};
            end
            
            if ~iscell(method)
                imethod = method;
                method  = cell(1,numel(signals));
                for is = 1:numel(signals)
                    method{is} = imethod;
                end
            end
            
            if ~exist('p0', 'var')
                p0=[];
            end
            
            this.signals = signals;
            this.params = {};
            this.method = method;
            this.p0 = p0;
            
            for ku = 1:numel(signals)
                this.params= [this.params { ...
                    [signals{ku} '_seed'],...
                    [signals{ku} '_dt_min'],...
                    [signals{ku} '_dt_max'],...
                    [signals{ku} '_min'],...
                    [signals{ku} '_max'],...
                    }];
                this.params_domain = [this.params_domain ... 
                                      BreachDomain('int') ...
                                      BreachDomain() ...
                                      BreachDomain() ...
                                      BreachDomain() ...
                                      BreachDomain() ...
                                      ];
            end
            
            if  isempty(this.p0)
                this.p0 = repmat( [0 0.1 1 -1 1] , 1, numel(signals));
            end
           
            % domains 
           this.signals_domain = repmat(BreachDomain(), 1, numel(this.signals));
            
        end
        
        function [X, time] = computeSignals(this,p, time)
            if size(p,1) ==1
                p = p';
            end
            
            X = zeros(numel(this.signals),numel(time));
            for isg = 1:numel(this.signals)
                
                rng(p((isg-1)*5+1), 'twister');    
                dt_min = p((isg-1)*5+2);
                dt_max = p((isg-1)*5+3);
                cp_min = p((isg-1)*5+4);
                cp_max = p(isg*5);
                
                t_cp = 0;
                while t_cp(end)<time(end)
                   dt_cp =  dt_min+(dt_max-dt_min)*rand(floor(10*time(end)/(dt_min+((dt_max-dt_min)/2))),1); 
                   t_cp = [t_cp;t_cp(end)+unique(cumsum(dt_cp))]; 
                end
                it_cp = find(t_cp>time(end)); 
                t_cp = t_cp(1:it_cp); 
                x_values = cp_min+ (cp_max-cp_min)*rand(1, numel(t_cp));
                
                meth = this.method{isg};
                
                if numel(t_cp)==1
                    x = x_values(end)*ones(numel(time),1);
                else
                    t_cp = t_cp(1:min([numel(t_cp) numel(x_values)]));
                    x = interp1(t_cp, x_values, time', meth, 'extrap');
                end
                X(isg,:) = x';
            end
        end
        
        function type = getType(this)
            type = 'varstep';
        end
        
        function args = getSignalGenArgs(this)
            args = {'method'};
        end
        
    end
end
