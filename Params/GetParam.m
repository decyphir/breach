function ParamValues = GetParam(S,ParamList)
% GETPARAM get the values of parameters in a parameter set
% 
% Synopsis: ParamValues = GetParam(P, ParamList)
%  
% Example (for Lorenz84 system):
%  
%    CreateSystem;
%    P = CreateSampling(Sys, {'a', 'b'}, [0 10; 0 5]);
%    Pr = Refine(P, 3);
%    val = GetParam(Pr, 'a' );
%    val = GetParam(Pr, 'b' );
%    val = GetParam(Pr, {'a', 'b'});
%    val = GetParam(Pr, {'b','a'});
%  
  
  
  if (isfield(S,'pts'))
    p= S.pts;
  else
    p = S.p;
  end
  
  if isstr(ParamList) 
    ind = FindParam(S,ParamList);    
    ParamValues = p(ind,:);
    return
  end
  
  if isnumeric(ParamList)
    ParamValues = [];  
    for i= 1:numel(ParamList)
      ParamValues(i,:) = p(ParamList(i),:);    
    end
    return
  end
    
  ParamValues = [];  
  for i= 1:numel(ParamList)
    ind = FindParam(S,ParamList{i});
    ParamValues(i,:) = p(ind,:);    
  end
    
  
