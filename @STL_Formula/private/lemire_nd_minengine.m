%*************************************************************************
% MATLAB MEX ROUTINE LEMIRE_ND_MINENGINE.C
% [minval idxmax] = LEMIRE_ND_MINENGINE(a, idx, window, shapeflag)
%
% PURPOSE: multiple 1D min running/filtering
% Similar to LEMIRE_ENGINE but working on the second dimension, while
% looping along the first and third dimensions. This MEX is used for
% the engine for multidimensional min/max filtering
%
% INPUTS
%  A: 3D arrays, logical and all numeric classes are supported
%  idx: 3D arrays, double, user input, must have the same number
%       of elements as A
%  window: scalar, size of the sliding window, must be >= 1
%  shapeflag: double scalar: 1, 2, 3 resp. for valid, same and full shape
% 
% OUTPUTS
%  For "valid" shape (without shapeflag passing)
%  minval: running min, vectors of dimension (length(A)-window+1), i.e.,
%      minval(:,1,:) is min(A(:,1:win,:))
%      minval(:,2,:) is min(A(:,2:win+1,:))
%      ...
%      minval(:,end,:) is min(A(:,end-window+1:end,:))
%  For "Same" shape output has the same dimension as A
%  For "Full" shape (with shapeflag)
%  output MINVAL has dimension (length(A)+window-1), correspond to all 
%  positions of sliding window that is intersect with A
%  
%  minidx: 3D arrays, subjected the same assignment as minval
%          The main purpose is to keep track of indexing during
%          runing filter
%
% Note: if the data is complex, the imaginary part is ignored.
%       window is limited to 2147483646 (2^31-2)
%
% Algorithm: Lemire's "STREAMING MAXIMUM-MINIMUM FILTER USING NO MORE THAN
% THREE COMPARISONS PER ELEMENT" Nordic Journal of Computing, Volume 13,
% Number 4, pages 328-339, 2006.
%
% Compilation:
%  >> mex -O -v lemire_nd_minengine.c 
% % add -largeArrayDims on 64-bit computer
%  >> mex -largeArrayDims -O -v lemire_nd_minengine.c
%
% see aldo: lemire_nd_maxengine.c, minmaxfilter
%           median filter, Kramer & Bruckner filter
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Contributor: Vaclav Potesil
% History
%  Original: 20/Sep/2009
%  Last update: 22/Sep/2009, input shapeflag always required
%                            same shape scan
%***********************************************************************/

error('Max file LEMIRE_ND_MINENGINE.C muts be compiled');
