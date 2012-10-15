function P = SetParam(P,ParamList,ParamValues)
% SETPARAM set the values of parameters in a parameter set
% 
% Synopsis: P = SetParam(P, ParamList, ParamValues)
% 
% Note that if the parameter is not present in P, it is created and
% appended.
% 
% Example (for Lorenz84 system):
% 
%    CreateSystem;
%    P = CreateSampling(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    val = GetParam(Pr, 'a');
%    val10 = 10*val; 
%    Pr10 = SetParam(Pr, 'a', val10); % values for 'a' in Pr10 are ten
%                                     % times those in Pr
% Other example :
% 
%    CreateSystem;
%    P = CreateSampling(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    Pr_2 = SetParam(Pr, {'F','G'}, [0.4; 0.6]); % values for 'F' (resp.
%                     % 'G') equals to 0.4 (resp 0.6) in the nine points.
%                     % Values of 'a' and 'b' are equals in Pr and Pr_2.
%    
%    
% SEE ALSO GETPARAM, CREATEPARAMSET 
%  

if size(ParamValues,1) == 1
    ParamValues = ParamValues';
end


if ischar(ParamList)
    ind = FindParam(P,ParamList);
    if (ind>size(P.pts,1))
        P.ParamList = [P.ParamList, {ParamList}];
    end
    P.pts(ind,:) = ParamValues;
    return

elseif iscell(ParamList)
    inds = FindParam(P,ParamList);
    new_params = ParamList(inds>size(P.pts,1));
    P.ParamList = [P.ParamList new_params];
    
    for ii= 1:numel(inds)
        P.pts(inds(ii),:) = ParamValues(ii,:);
    end
    return

elseif isnumeric(ParamList)
    for i = 1:numel(ParamList)
        if ParamList(i) > size(P.pts, 1)
            error('SetParam:index','Index out of range')
        end
        P.pts(ParamList(i),:) = ParamValues(i,:);
    end
    return
end

% Why the following code ??
for i = 1:numel(ParamList)
    ind = FindParam(P,ParamList{i});
    P.pts(ind,:) =  ParamValues(i,:);
end

end