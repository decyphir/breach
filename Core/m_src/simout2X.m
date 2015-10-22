function [tout, X, signals] = simout2X(simout)
%
% converts a simulink output to a data structure Breach can handle
%

tout = simout.get('tout')';
X=[];
%% Outputs and scopes
Vars = simout.who;
lenVars = numel(Vars);
signals = {};

for iV = 1:lenVars
    Y = get(simout,Vars{iV});
    if ~isempty(Y)
        
        if ~strcmp(Vars{iV}, 'tout')&&~strcmp(Vars{iV},'logsout')&&(isstruct(Y))
            for iS=1:numel(Y.signals)
                signame = Y.signals(iS).label;
                if ~ismember(signame,signals)
                    
                    nbdim = size(double(Y.signals(iS).values),2);
                    try
                        xx = interp1(Y.time, double(Y.signals(iS).values),tout, 'linear','extrap') ;
                    catch
                        if (nbdim==1)
                            xx = 0*tout;
                        else
                            xx = zeros(numel(tout), nbdim);
                        end
                    end
                    
                    if (nbdim==1)
                        X = [X; xx];
                        signals = {signals{:} signame };
                    else
                        X = [X; xx'];
                        for idim = 1:nbdim
                            signamei = [signame '_' num2str(idim)  '_'];
                            signals = {signals{:} signamei};
                        end
                    end
                end
            end
        end
    end
end

logs = simout.get('logsout');

if ~isempty(logs)
    logs_names = logs.getElementNames();
    
    %% logs
    for ilg = 1:numel(logs_names)
        if ~(ismember(logs_names{ilg}, signals))
            signame = logs_names{ilg};
            if ~ismember(signame,signals)
                
                sig = logs.getElement(signame);
                nbdim = size(sig.Values.Data,2);
                
                % getting signal data
                for idim =1:nbdim
                    try
                        xdata = interp1(sig.Values.Time',double(sig.Values.Data(:,idim)),tout, 'linear','extrap');
                        X = [X ; xdata(1,:)];
                    end
                end
                
                % naming multidimensional signal= name_signal_i_
                if nbdim==1
                    signals = {signals{:} signame};
                else
                    for idim =1:nbdim
                        signamei = [signame '_' num2str(idim)  '_'];
                        signals = {signals{:} signamei};
                    end
                end
            end
        end
    end
end
end
