function [varargout] = minmaxfilt1(A, window, outtype, shape)
% SHORT CALL: [minval, maxval] = MINMAXFILT1(A, window)
% 
% LONG CALL: [minval, maxval] = MINMAXFILT1(A, window, outtype, shape)
%
% PURPOSE: one-dimensional min/max (running) filtering
%
%  INPUTS
%   A: 1D array, logical and all numeric classes are supported
%   WINDOW: size of the sliding window. The value must be >= 1
%           - scalar -> same window size will be scanned for all dimension
%   OUTTYPE: ['both'], 'minmax', 'min', 'max'. If outtype is 'both' or
%            'minmax', both MIN and MAX arrays will be returned (in that
%            order). Otherwise MINMAXFILT1 returns only the requested array.
%   SHAPE: 'full' 'same' ['valid']
%       'full'  - (default) returns the full size arrays,
%       'same'  - returns the central part of the filtering output
%                 that is the same size as A.
%       'valid' - returns only those parts of the filter that are computed
%                 without the padding edges.
%  OUTPUTS
%   minval, maxval: filtered min/max arrays of same number of 
%   dimension as the input A. The size depends on SHAPE. If SHAPE is
%       - 'full', the size is size(A)+(WINDOW-1)
%       - 'same', the size is size(A)
%       - 'valid', the size is size(A)-WINDOW+1
%
% [MINVAL MAXVAL MINIDX MAXIDX] = MINMAXFILT1(A, ...) returns index arrays
%   such that:
%       A(MINIDX) is equal to MINVAL and
%       A(MAXIDX) is equal to MAXVAL
%   
%  Note: if the data is complex, the imaginary part is ignored.
%
%  Algorithm: Lemire's "STREAMING MAXIMUM-MINIMUM FILTER USING NO MORE THAN
%  THREE COMPARISONS PER ELEMENT" Nordic Journal of Computing, Volume 13,
%  Number 4, pages 328-339, 2006.
%
% See also: minmaxfilt, ordfilt2, medfilt2 (Image Processing Toolbox)
%
% AUTHOR: Bruno Luong <brunoluong@yahoo.com>
% Contributor: Vaclav Potesil
% HISTORY
%   Original: 12-Jul-2009
%   20/Sep/2009, optional output runing min/max indexing
%                 use ND engine (to ease code maintenance)
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
    error('minmaxfilt1:badOuttype','MINMAXFILT1: unknown outtype %s', outtype);
end
%NM: suggestion to avoid use of strmatch which is obsolete:
%try
%    validatestring(outtype,{'both' 'minmax' 'maxmin' 'min' 'max'});
%catch %#ok<CTCH>
%    error('minmaxfilt1:badOuttype','MINMAXFILT1: unknown outtype %s', outtype);
%end

% Default output types, both min and max
if nargin<4 || isempty(shape)
    shape = 'valid'; % shape
end
shape = lower(strtrim(shape));
% Check if SHAPE is correct
shapeloc = strmatch(shape,{'valid' 'same' 'full'});
if isempty(shapeloc)
    error('minmaxfilt1:unknownShape','MINMAXFILT1: unknown shape %s', shape);
end
%NM: suggestion to avoid use of strmatch which is obsolete:
%try
%    validStr = {'valid' 'same' 'full'};
%    shape = validatestring(shape,validStr);
%    shapeloc = find(shape,validStr);
%catch %#ok<CTCH>
%    error('minmaxfilt1:unknownShape','MINMAXFILT1: unknown shape %s', shape);
%end

window = min(window,numel(A)-1);

isrow = size(A,1)==1;
% reshape in row
A = reshape(A, 1, []);
n = size(A,2);

% Output
out = cell(1,0);
idx = cell(1,0);
nout = 0;

% MIN
if ~strcmp(outtype,'max')
    % Index array, they must be double
    minidx = 1:double(n);
    % Call MEX engines
    [minval, minidx] = lemire_nd_minengine(A, minidx, window, shapeloc);
    % reshape in row
    if ~isrow
        minval = minval(:);
        minidx = minidx(:);
    end
    nout=nout+1;
    out{nout} = minval;
    idx{nout} = minidx;
end

% MAX
if ~strcmp(outtype,'min')
    % Index array, they must be double
    maxidx = 1:double(n);
    % Call MEX engines
    [maxval, maxidx] = lemire_nd_maxengine(A, maxidx, window, shapeloc);
    % reshape in row
    if ~isrow
        maxval = maxval(:);
        maxidx = maxidx(:);
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

end % minmaxfilt1
