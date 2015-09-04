function [varargout] = minmaxfilt(A, window, outtype, shape)
% SHORT CALL: [minimg maximg] = MINMAXFILT(A, window)
% 
% LONG CALL: [minimg maximg] = MINMAXFILT(A, window, outtype, shape)
%
% PURPOSE: multi-dimensional min/max (erosion/dilatation) filtering
%
%  INPUTS:
%   A: ND arrays, logical and all numeric classes are supported
%   WINDOW: size of the sliding window. The value must be >= 1
%           - scalar -> same window size will be scanned for all dimension
%           - [n1, ...,nd] -> define separate scan size for each dimension
%   OUTTYPE: ['both'], 'minmax', 'min', 'max'. If outtype is 'both' or
%            'minmax', both MIN and MAX arrays will be returned (in that
%            order). Otherwise MINMAXFILT returns only the requested array.
%   SHAPE: 'full' 'same' ['valid']
%       'full'  - (default) returns the full size arrays,
%       'same'  - returns the central part of the filtering output
%                 that is the same size as A.
%       'valid' - returns only those parts of the filter that are computed
%                 without the padding edges.
%  OUTPUTS:
%   minimg, maximg: filtered min/max arrays of same number of 
%   dimension as the input A. The size depends on SHAPE. If SHAPE is
%       - 'full', the size is size(A)+(WINDOW-1)
%       - 'same', the size is size(A)
%       - 'valid', the size is size(A)-WINDOW+1
%
% [MINIMG MAXIMG MINIDX MAXIDX] = MINMAXFILT(A, ...) returns the linear
% index arrays such that:
%       A(MINIDX) is equal to MINIMG and
%       A(MAXIDX) is equal to MAXIMG
%   
%  Note: if the data is complex, the imaginary part is ignored.
%
%  Algorithm: Lemire's "STREAMING MAXIMUM-MINIMUM FILTER USING NO MORE THAN
%  THREE COMPARISONS PER ELEMENT" Nordic Journal of Computing, Volume 13,
%  Number 4, pages 328-339, 2006.
%
% See also: minmaxfilt1, ordfilt2, medfilt2 (Image Processing Toolbox)
%
% AUTHOR: Bruno Luong <brunoluong@yahoo.com>
% Contributor: Vaclav Potesil
% HISTORY
%   Original: 12-Jul-2009
%   Last update:
%   20/Sep/2009, optional output runing min/max indexing
%                separated min/max and improved engines
%                correct a bug of crop indexes (reported by 
%                Vaclav Potesil)
%   22/Sep/2009, simplify code when handling same shape

% Default sliding window size
if nargin<2 || isempty(window)
    window = 3;
end

% Default output types, both min and max
if nargin<3 || isempty(outtype)
    outtype = 'both'; % min/max
end
outtype = lower(strtrim(outtype));
% Check if OUTTYPE is correct
if isempty(strmatch(outtype,{'both' 'minmax' 'maxmin' 'min' 'max'}))
    error('MINMAXFILT: unknown outtype %s', outtype);
end

% Default output types, both min and max
if nargin<4 || isempty(shape)
    shape = 'valid'; % shape
end
shape = lower(strtrim(shape));
% Check if SHAPE is correct
shapeloc = strmatch(shape,{'valid' 'same' 'full'});
if isempty(shapeloc)
    error('MINMAXFILT: unknown shape %s', shape);
end

% We do not support SPARSE
if issparse(A)
    error('MINMAXFILT: first output A must be full matrix')
end

nd = ndims(A);
% extend window size
if isscalar(window) % same for all dimensions
    window(1:nd) = window;
else
    % pad missing window size
    window(end+1:nd) = 1;
end

out = cell(1,0);
idx = cell(1,0);
nout = 0;

szA = size(A);
    
% MINVAL
if ~strcmp(outtype,'max')
    minval = A;
    % Create linear index array
    minidx = zeros(szA,'double');
    minidx(:) = 1:numel(A);    
    sz = size(minval);
    % Loop on dimension
    for dim=1:nd
        % Reshape minval in 3D arrays, working dimension is the middle
        % That is the form required by LEMIRE_ND_ENGINE
        p = prod(sz(1:dim-1)); % return 1 if empty
        n = sz(dim);
        q = prod(sz(dim+1:end)); % return 1 if empty
        minval = reshape(minval,[p n q]);
        win = window(dim);
        % call mex engine
        if win~=1
            [minval minidx] = lemire_nd_minengine(minval, minidx, win, shapeloc);
        end
        % Blow back to n-dimensional
        sz(dim) = size(minval,2);
        minval = reshape(minval, sz);
        minidx = reshape(minidx, sz);
    end
    nout=nout+1;
    out{nout} = minval;
    idx{nout} = minidx;
end

% MAXVAL
if ~strcmp(outtype,'min')
    maxval = A;
    % Create linear index array
    maxidx = zeros(szA,'double');
    maxidx(:) = 1:numel(A);
    sz = size(maxval);
  
    % Loop on dimension
    for dim=1:nd
        % Reshape maxval in 3D arrays, working dimension is the middle
        % That is the form required by LEMIRE_ND_ENGINE
        p = prod(sz(1:dim-1)); % return 1 if empty
        n = sz(dim);
        q = prod(sz(dim+1:end)); % return 1 if empty
        maxval = reshape(maxval,[p n q]);
        win = window(dim);
        % call mex engine
        if win~=1
            [maxval maxidx] = lemire_nd_maxengine(maxval, maxidx, win, shapeloc);
        end
        % Blow back to n-dimensional
        sz(dim) = size(maxval,2);       
        maxval = reshape(maxval, sz);
        maxidx = reshape(maxidx, sz);
    end
    nout=nout+1;
    out{nout} = maxval;
    idx{nout} = maxidx;
end

% Assign the output
if nargout <= length(out)
    varargout = out;
else
    out = [out idx];
    varargout = out;
end

end % minmaxfilt
