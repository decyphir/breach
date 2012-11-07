function P = SAddUncertainParam(P, is)
% SADDUNCERTAINPARAM add uncertain parameters to a parameter set. The epsi
% for these new parameters is set to 1/10th of the value of the parameter
% if it is not null, to zero otherwize.
%
%  Synopsis : P = SAddUncertainParam(P, is)
%
%  Input :
%    - P  The parameter set to change
%    - is The index of the parameter(s) to add
%
%  Output :
%     - P  The modified parameter set. In this parameter set, all the
%            parameters indexed by is are considered as uncertain.
%
%  Examples for Lorenz84 :
%    CreateSystem;
%    P = CreateParamSet(Sys,{'a'},[0.1 0.8]);
%    P = SAddUncertainParam(P, FindParam(P,{'b','F'}));
%

if(~isempty(find(is<0,1)))
   return ; 
end

if(~isempty(find(is>numel(P.ParamList),1)))
    return ;
end

epsi = zeros(1,size(P.pts,2));
for i = is
    if (isempty(P.dim(P.dim==i)))
        P.dim = [P.dim i];
        ptsun = P.pts(i,:); %on récupère les valeurs du paramètre
        epsi(ptsun~=0) = abs(ptsun(ptsun~=0)/10);
        epsi(ptsun==0) = 1;
        P.epsi(end+1,:) = epsi;
    end
end