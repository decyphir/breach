function S = SetParam(S,ParamList,ParamValues)
% SETPARAM set the values of parameters in a parameter set
% 
% Synopsis: ParamValues = SetParam(P, ParamList, ParamValues)
%  
% Note that if the parameter is not present in P, it is created and appended.  
%  
% Example (for Lorenz84 system):
%  
%    CreateSystem;
%    P = CreateSampling(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    val = GetParam(Pr, 'a' );
%    val10 = 10*val; 
%    Pr10 = SetParam(Pr, 'a', val10); % values for 'a' in Pr10 are ten times those in Pr
%  
% SEE ALSO GETPARAM, CREATEPARAMSET 
%  
  
  if isstr(ParamList) 
    ind = FindParam(S,ParamList);
    S.ParamList = unique({S.ParamList{:}, ParamList},'stable');
    S.pts(ind,:) = ParamValues;
    return    
  elseif iscell(ParamList)
    ind = FindParam(S,ParamList);
    S.ParamList = unique({S.ParamList{:}, ParamList{:}},'stable');
    S.pts(ind,:) = ParamValues;
    return        
  end
  
  if isnumeric(ParamList)
    for i= 1:numel(ParamList)
      S.pts(ParamList(i),:) = ParamValues(i,:);
    end
    return
  end
  
  for i= 1:numel(ParamList)
    ind = FindParam(S,ParamList{i});
    S.pts(ind,:) =  ParamValues(i,:);    
  end
    
