function P = SAddUncertainParam(P, param)
% SADDUNCERTAINPARAM add uncertain parameters to a parameter set. The epsi
% for these new parameters is set to 1/10th of the value of the parameter
% if it is not zero, to one otherwise.
%
%  Synopsis : P = SAddUncertainParam(P, param)
%
%  Input :
%    - P     : The parameter set to modify
%    - param : The indexes or names of the parameter(s) to add. If any
%              param is not valid (ie: negative index or index higher than
%              size(P.pts,1) or param not in ParamList), nothing is done.
%
%  Output :
%     - P  The modified parameter set. In this parameter set, all the
%          parameters indexed by idx are considered as uncertain.
%
%  Examples (Lorenz84):
%    CreateSystem;
%    P = CreateParamSet(Sys,{'a'},[0.1 0.8]);
%    
%    P1 = SAddUncertainParam(P, {'b','F'});
%    P1.ParamList(P1.dim)  % shows uncertain parameters
%
%    P2 = SAddUncertainParam(P, 'blah');  % 'blah' does not exist
%    P2.ParamList(P2.dim)  % shows uncertain parameters
%    P2 = SAddUncertainParam(P, 'x0');
%    P2.ParamList(P2.dim)  % shows uncertain parameters
%
%See also SDelUncertainParam CreateParamSet SetEpsi GetEpsi
%

if(ischar(param)||iscell(param))
    param = FindParam(P,param);
end

if(~isempty(find(param<0,1)) || ~isempty(find(param>size(P.pts,1),1)))
    warning('SAddUncertainParam:NotValidParam','A parameter name or index is not correct');
   return ; 
end

epsi = zeros(1,size(P.pts,2));
for ii = param
    if ~ismember(ii,P.dim)
        P.dim = [P.dim ii];
        ptsun = P.pts(ii,:);  % we gather the parameter values
        epsi(ptsun~=0) = abs(ptsun(ptsun~=0)/10);
        epsi(ptsun==0) = 1;
        P.epsi(end+1,:) = epsi;
    end
end

end
