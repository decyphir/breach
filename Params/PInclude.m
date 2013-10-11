function [include, info] = PInclude(Pin, P, EPSI)
%PINCLUDE tests if the parameter set Pin is included in P. Pin.ParamList
% must be included in P.ParamList. This function tests if all parameter
% vectors in Pin are included in, at least, one parameter vector of P. If a
% parameter pp is a fixed parameter in P and in Pin, Pin is said included
% in P iff
%   GetParam(P,'pp')-EPSI <= GetParam(Pin,'pp') <= GetParam(P,'pp')+EPSI
% If the parameter pp is an uncertain parameter of Pin, Pin is said
% included in P iff
%   GetParam(P,'pp')-EPSI <= GetParam(Pin,'pp')-GetEpsi(Pin,'pp)  &&
%   GetParam(Pin,'pp')+GetEpsi(Pin,'pp') <= GetParam(P,'pp')+EPSI
%
% Synopsis: [include, info] = PInclude(Pin, P [, EPSI] )
%
% Input:
%  - Pin  : A parameter set which may contain many parameter vectors
%  - P    : A parameter set with the same parameters than Pin
%  - EPSI : optional, default=eps
%
% Output:
%  - include : true if all parameter sets in Pin are included in a
%              parameter set of P, false otherwise.
%  - info    : provides more details if Pin is not included in P. If
%              Pin.ParamList is not included in P.ParamList, info is set to
%              -1. Otherwise, info indicates the indexes of the parameter
%              vectors in Pin which are not included in P. If Pin is
%              included in P, info is set to 0.
%
%See also SConcat SSelect
%

include = true;

if ~isempty(setdiff(Pin.ParamList,P.ParamList)) % we check that all parameter in Pin are in P
    include = false;
    info = -1;
    return;
end

numPparam = numel(P.ParamList); % number of parameters
if(numPparam~=numel(Pin.ParamList))
    include = false;
    info = -1;
    return;
end

if(nargin<=2)
    EPSI = eps;
end

numPdim = numel(P.dim);

% we create a 3D array describing the intervals in the following order:
% uncertains param of P in the order defined by P.dim, then,
% fixed parameters of P in the order defined by P.ParamList
Pinterv = zeros([size(P.pts),2]);
Pinterv(1:numPdim,:,1) = P.pts(P.dim,:) - P.epsi; % first component in the 3rd dim = min value
Pinterv(1:numPdim,:,2) = P.pts(P.dim,:) + P.epsi; % second component in the 3rd dim = max value
Pinterv(numPdim+1:end,:,1) = P.pts(setdiff(1:numPparam,P.dim,'stable'),:);
Pinterv(numPdim+1:end,:,2) = Pinterv(numPdim+1:end,:,1); % fixed param -> min=max

% we compute the order for Pin
paramOrder = FindParam(Pin,P.ParamList([P.dim,setdiff(1:numPparam,P.dim,'stable')]));
% paramOrder is such that Pin.pts(paramOrder,:) is in the same order than
% Pinterv
[~,uncertainParamOrder] = intersect(paramOrder,Pin.dim,'stable');
% uncertainParamOrder is such that Pinterv(uncertainParamOrder,:,:) correspond
% to the same param than Pin.pts(Pin.dim,:)
[fixedParamOrderPts,fixedParamOrderVal] = setdiff(paramOrder,Pin.dim,'stable');
% fixedParamOrderPts and fixedParamOrderVal are such that
% Pin.pts(fixedParamOrderPts,:) correspond to the same param than
% Pinterv(fixedParamOrderVal,:,:)

info = [];
for ii = 1:size(Pin.pts,2)
    PinValue = zeros(numPparam,1,2); % PinValue contains the intervals for the ii-th param set in Pin
    PinValue(uncertainParamOrder,1,1) = Pin.pts(Pin.dim,ii) - Pin.epsi(:,ii);
    PinValue(uncertainParamOrder,1,2) = Pin.pts(Pin.dim,ii) + Pin.epsi(:,ii);
    PinValue(fixedParamOrderVal,1,1) = Pin.pts(fixedParamOrderPts,ii);
    PinValue(fixedParamOrderVal,1,2) = PinValue(fixedParamOrderVal,1,1);
    
    PinValue = repmat(PinValue,[1,size(P.pts,2),1]); % replicate PinValue, so it has the same size than Pinterv
    
    if ~any(all(Pinterv(:,:,1)-EPSI <= PinValue(:,:,1) & PinValue(:,:,2) <= Pinterv(:,:,2)+EPSI , 1))
        info = [info,ii]; %#ok<AGROW>
        include = false;
%        return;
    end
    
end

if(include)
    info = 0;
end

end

