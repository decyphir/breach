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

for ii = 1:numel(P.ParamList)
    param = P.ParamList{ii};
    val = GetParam(P,param);
    spacing = repmat(' ',1,12-numel(param));
    mini = min(val);
    maxi = max(val);
    if(abs(maxi-mini)<eps)
        fprintf(fid,'%s%s :  %.2g\n',param,spacing,mini);
    else
        fprintf(fid,'%s%s :  %.2g -- %.2g\n',param,spacing,mini,maxi);
    end
end

end