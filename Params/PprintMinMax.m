function PprintMinMax(P, file, access)
%PPRINTMINMAX writes the extrema of the parameter vectors of P.
%
% Synopsis: PprintMinMax(P[, file[, access]])
% Inputs:
%  - P      : the parameter set. It may contain many parameter vector
%  - file   : (optional, default print to the screen) indicates the file in
%             which the extrema should be written.
%  - access : (optional, default = 'w') indicates if the output must erase
%             the file (in which case, access must be equal to 'w') or be
%             appended at the end of the file (in which case, access must
%             be equal to 'a').
%
% Output:
%  - None
%

if(~exist('file','var') || isempty(file))
    fid = 1; % standart output (ie: the screen)
else
    if(~exist('access','var') || (access~='w' && access~='a'))
        access = 'w';
    end
    [fid,err_msg] = fopen(file,access);
    if(fid==-1)
        warning('PprintMinMax:fileAccesFailed',['The file ''',file,'''couldn''t be opened. Error message is:\n',err_msg] );
    end
end
fprintf(fid,'\nStatistics based on %d parameter vectors.\n\n',size(P.pts,2));
fprintf(fid,'\nname         :  interval            |');


% Etude du nombre de centres
fprintf(fid,'nb center\n');
figure;

% % Analyse des données selon une loi log-normale
% fprintf(fid,'mean        std. dev. \n');


fprintf(fid,[repmat('-',1,61) '\n']);


for ii = 1:numel(P.ParamList)
    if(mod(ii,3)==1 && ii~=1)
        fprintf(fid,[repmat('-',1,61) '\n']);
    end
    param = P.ParamList{ii};
    val = GetParam(P,param);
    mini = min(val);
    maxi = max(val);
    
    
    % Etude du nombre de centres
    val_log = log(val);
    val_log(val_log==-Inf) = realmin;
    [~,C_old,sumd_old] = kmeans(val_log,1);
    for jj=2:10
        [~,C,sumd] = kmeans(val_log,jj,'emptyaction','drop');
        if(sum(sumd_old)/sum(sumd)<jj/(jj-1)) % jj/(jj-1) correspond au cas homogène ; on ajoute un facteur 2 parce que c'est cool
            %sumd = sumd_old;
            C = sort(C_old);
            break
        else
            sumd_old = sumd;
            C_old = C;
        end
    end
    subplot(numel(P.ParamList),1,ii);
    plot(val,ones(1,numel(val)),'*');
    set(gca,'XScale','log','YTick',[]);
    ylabel(P.ParamList{ii},'Interpreter','none','Rotation',0,'HorizontalAlignment','Right');
    
%     % Analyse des données selon une loi log-normale
%     if(abs(maxi)<=eps && abs(mini)<eps)
%         M = 0;
%         V = 0;
%         %parmci = [0,0;0,0];
%     else
%         %[parmhat,parmci] = lognfit(val);
%         parmhat = lognfit(val);
%         [M,V] = lognstat(parmhat(1),parmhat(2));
%     end
    

    spacing_param = repmat(' ',1,12-numel(param));
    if(abs(maxi-mini)<eps)
        val_str = sprintf('%.2g',mini);
    else
        val_str = sprintf('%.2g -- %.2g',mini,maxi);
    end
    spacing_val = repmat(' ',1,18-numel(val_str));
    fprintf(fid,'%s%s :  %s%s  | ',param,spacing_param,val_str,spacing_val);
    
    % Etude du nombre de centres
    fprintf(fid,'%d center: ',numel(C));
    for kk=1:numel(C)
        fprintf(fid,'%.2g ',exp(C(kk)));
    end
    fprintf(fid,'\n');
    
%     % Analyse des données selon une loi log-normale
%     mean_str = sprintf('%.2g',M);
%     spacing_mean = repmat(' ',1,8-numel(mean_str));
%     fprintf(fid,'E: %s%s SD: %.2g\n',mean_str,spacing_mean,sqrt(V))
end

end
