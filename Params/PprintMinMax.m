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
    fid = fopen('file',access);
end

fprintf(fid,'\n nom         :  valeur              | moyenne     ecart type\n');
fprintf(fid,[repmat('-',1,61) '\n']);

for ii = 1:numel(P.ParamList)
    if(mod(ii,3)==1 && ii~=1)
        fprintf(fid,[repmat('-',1,61) '\n']);
    end
    param = P.ParamList{ii};
    val = GetParam(P,param);
    mini = min(val);
    maxi = max(val);
    if(abs(maxi)<=eps && abs(mini)<eps)
        M = 0;
        V = 0;
        %parmci = [0,0;0,0];
    else
        %[parmhat,parmci] = lognfit(val);
        parmhat = lognfit(val);
        [M,V] = lognstat(parmhat(1),parmhat(2));
    end
    
    spacing_param = repmat(' ',1,12-numel(param));
    if(abs(maxi-mini)<eps)
        val_str = sprintf('%.2g',mini);
    else
        val_str = sprintf('%.2g -- %.2g',mini,maxi);
    end
    spacing_val = repmat(' ',1,18-numel(val_str));
    mean_str = sprintf('%.2g',M);
    spacing_mean = repmat(' ',1,8-numel(mean_str));
    fprintf(fid,'%s%s :  %s%s  | E: %s%s SD: %.2g\n',param,spacing_param,val_str,spacing_val,mean_str,spacing_mean,sqrt(V));
end

end