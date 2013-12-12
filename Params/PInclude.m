function [include, outside] = PInclude(Pin, P, epsi)
%PINCLUDE tests if the parameter set Pin is included in P. Pin and P must
% contains the same parameters, but they may be in different order. This
% function tests if all parameter vectors in Pin are included in, at least,
% one parameter vector of P. If a parameter pp is a fixed parameter in P
% and in Pin, then Pin is said included in P iff
%   GetParam(P,'pp')*(1-epsi) <= GetParam(Pin,'pp') <= GetParam(P,'pp')*(1+epsi)
% If the parameter pp is an uncertain parameter of Pin and fixed in P, Pin
% is said included in P iff
%   GetParam(P,'pp')*(1-epsi) <= GetParam(Pin,'pp')-GetEpsi(Pin,'pp')  &&
%   GetParam(Pin,'pp')+GetEpsi(Pin,'pp') <= GetParam(P,'pp')*(1+epsi)
% 
% Synopsis: [include, outside] = PInclude(Pin, P [, epsi] )
% 
% Inputs:
%  - Pin  : A parameter set which may contain many parameter vectors
%  - P    : A parameter set with the same parameters than Pin
%  - epsi : Optional, default=eps. Indicate the admissible relative error.
% 
% Outputs:
%  - include : true if all parameter sets in Pin are included in a
%              parameter set of P, false otherwise.
%  - info    : provides more details if Pin is not included in P. If
%              Pin.ParamList is not included in P.ParamList, info is set to
%              -1. Otherwise, info indicates the indexes of the parameter
%              vectors in Pin which are not included in P. If Pin is
%              included in P, info is set to 0.
% 
% Example (Lorentz84):
%   CreateSystem
%   P = CreateParamSet(Sys,{'x0','x1','a'},[-5,-2;3,4;0,0.5]);
%   Pin = CreateParamSet(Sys,{'x1','x0'},[3.5,3.8;-4,-3]);
%   PInclude(Pin,P) % should be 1
%   
%   Pout = SetParam(Pin,'a',1);
%   PInclude(Pout,P) % should be 0
%   
%   Pin = SAddUncertainParam(Pin,'a');
%   Pin = SetParam(Pin,'a',0.35);
%   Pin = SDelUncertainParam(Pin,'x0');
%   Pin = SAddUncertainParam(Pin,'x0'); % change order in epsi
%   PInclude(Pin,P)
%   
%   Pout = SetEpsi(Pin,'a',0.2);
%   PInclude(Pout,P)
%   
%   Pboth = Refine(Pout,3); % 27 parameter vectors (param 'a' is splitted
%                           % into [0.15, 0.283], [0.283, 0.416] and
%                           % [0.416 0.55] )
%   [include,outside] = PInclude(Pboth,P)
%   a = GetParam(Pboth,'a')+GetEpsi(Pboth,'a'); % upper bound for 'a'
%   a(outside) % should only get 0.55
% 
%See also SConcat SSelect SAddUncertainParam
%

include = true;

if ~all(ismember(Pin.ParamList,P.ParamList)) % we check that all parameter in Pin are in P
    include = false;
    outside = -1;
    return;
end

numPparam = numel(P.ParamList); % number of parameters
if(numPparam~=numel(Pin.ParamList)) % check that number of parameter is the same in P and Pin
    include = false;
    outside = -1;
    return;
end

if(nargin<=2)
    epsi = eps;
end

numPdim = numel(P.dim);

% we create a 3D array describing the intervals in the following order:
% uncertains param of P in the order defined by P.dim, then,
% fixed parameters of P in the order defined by P.ParamList
rangeP = zeros([size(P.pts),2]);
rangeP(1:numPdim,:,1) = P.pts(P.dim,:)*(1-epsi) - P.epsi; % first component in the 3rd dim = min value
rangeP(1:numPdim,:,2) = P.pts(P.dim,:)*(1+epsi) + P.epsi; % second component in the 3rd dim = max value
idxPParamFixed = setdiff(1:numPparam,P.dim,'stable'); % indexes of fixed parameter (ie not uncertain)
rangeP(numPdim+1:end,:,1) = P.pts(idxPParamFixed,:)*(1-epsi);
rangeP(numPdim+1:end,:,2) = P.pts(idxPParamFixed,:)*(1+epsi); % fixed param, so min=max

[~,uncertainParamOrder] = ismember(Pin.ParamList(Pin.dim),P.ParamList([P.dim,idxPParamFixed]));
uncertainParamOrder = uncertainParamOrder(uncertainParamOrder~=0);
% uncertainParamOrder is such that Pin.pts(Pin.dim,:) correspond to the
% same param than rangeP(uncertainParamOrder,:,:)
idxPinParamFixed = setdiff(1:numPparam,Pin.dim,'stable'); % indexes of fixed parameter (ie not uncertain)
[~,fixedParamOrder] = ismember(Pin.ParamList(idxPinParamFixed),P.ParamList([P.dim,idxPParamFixed]));
fixedParamOrder = fixedParamOrder(fixedParamOrder~=0);
% fixedParamOrder is such that Pin.pts(idxPinParamFixed,:) correspond to
% the same param than rangeP(fixedParamOrder,:,:)

outside = [];
for ii = 1:size(Pin.pts,2)
    rangePin = zeros(numPparam,1,2); % rangePin contains the intervals for the ii-th param set in Pin
    rangePin(uncertainParamOrder,1,1) = Pin.pts(Pin.dim,ii) - Pin.epsi(:,ii);
    rangePin(uncertainParamOrder,1,2) = Pin.pts(Pin.dim,ii) + Pin.epsi(:,ii);
    rangePin(fixedParamOrder,1,1) = Pin.pts(idxPinParamFixed,ii);
    rangePin(fixedParamOrder,1,2) = rangePin(fixedParamOrder,1,1); % fixed param, so min=max
    
    rangePin = repmat(rangePin,[1,size(P.pts,2),1]); % replicate rangePin, so it has the same size than rangeP
    
    if ~any(all(rangeP(:,:,1) <= rangePin(:,:,1) & rangePin(:,:,2) <= rangeP(:,:,2) , 1))
        outside = [outside,ii]; %#ok<AGROW>
        include = false;
    end
    
end

if(include)
    outside = 0;
end

end
