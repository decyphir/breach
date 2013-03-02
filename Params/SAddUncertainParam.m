function P = SAddUncertainParam(P, idx)
% SADDUNCERTAINPARAM add uncertain parameters to a parameter set. The epsi
% for these new parameters is set to 1/10th of the value of the parameter
% if it is not zero, to one otherwise.
%
%  Synopsis : P = SAddUncertainParam(P, idx)
%
%  Input :
%    - P  : The parameter set to change
%    - idx: The index of the parameter(s) to add. If any index in idx is
%           negative or if any index does not correspond to an existing
%           parameter, nothing is done.
%
%  Output :
%     - P  The modified parameter set. In this parameter set, all the
%          parameters indexed by idx are considered as uncertain.
%
%  Examples (Lorenz84):
%    CreateSystem;
%    P = CreateParamSet(Sys,{'a'},[0.1 0.8]);
%    P = SAddUncertainParam(P, FindParam(P,{'b','F'}));
%
%See also SDelUncertainParam
%

if ~isempty(find(idx<0,1))
   return ; 
end

if ~isempty(find(idx>size(P.pts,1),1))
    return ;
end

epsi = zeros(1,size(P.pts,2));
for ii = idx
    if ~ismember(ii,P.dim)
        P.dim = [P.dim ii];
        ptsun = P.pts(ii,:);  % we gather the parameter values
        epsi(ptsun~=0) = abs(ptsun(ptsun~=0)/10);
        epsi(ptsun==0) = 1;
        P.epsi(end+1,:) = epsi;
    end
end

end
