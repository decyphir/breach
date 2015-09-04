%  [minval maxval] = LEMIRE_ENGINE(a, window)
%  
%  PURPOSE: single shoot (one vector) 1D min/max running/filtering
%  
%  INPUTS
%   A: vector, logical and all numeric classes are supported
%   window: scalar, size of the sliding window, must be >= 1
%  
%  OUTPUTS
%   minval, maxval: running min/max, vectors of dimension
%   (length(A)-window+1), i.e.,
%       minval(1) is min(A(1:win))
%       minval(2) is min(A(2:win+1))
%       ...
%       minval(end) is min(A(end-window+1:end))
%   The same indexing arrangement applies for maxval
%  
%  Note: if the data is complex, the imaginary part is ignored.
%        This function has less overhead than LEMIRE_ND_ENGINE
%  
%  Algorithm: Lemire's "STREAMING MAXIMUM-MINIMUM FILTER USING NO MORE THAN
%  THREE COMPARISONS PER ELEMENT" Nordic Journal of Computing, Volume 13,
%  Number 4, pages 328-339, 2006.
%  
%  Compilation:
%   >> mex -O -v lemire_engine.c % add -largeArrayDims on 64-bit computer
%  
%  see aldo: lemire_idxengine, median filter, Kramer & Bruckner filter
%  
%  Author: Bruno Luong <brunoluong@yahoo.com>
%  History
%   Original: 12/July/2009 
