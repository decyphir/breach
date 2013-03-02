function P = SDelUncertainParam(P, idx, dims)
%SDELUNCERTAINPARAM removes uncertains parameters indexed by is in P
% excepted of those contained in dims.
%
%Synopsis : P = SDelUncertainParam(P, is, dims)
%
%See also SAddUncertainParam
%

% does not remove dims

if ~exist('dims','var')
    dims = [];
end

idx_remove = setdiff(idx,dims,'stable'); % index of uncertains parameter to remove
idx = setdiff(P.dim,idx_remove,'stable'); % index of uncertains parameter to keep
P.dim = P.dim(idx);
P.epsi = P.epsi(idx,:);

% for ii=idx
%     if isempty(find(dims==ii,1))
%         idx = find(P.dim~=ii);
%         P.dim = P.dim(idx);
%         P.epsi = P.epsi(idx,:);
%     end
% end

end
