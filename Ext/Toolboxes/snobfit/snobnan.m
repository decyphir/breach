% Copyright (c) 2003-2008, Arnold Neumaier
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the University of Vienna nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY ARNOLD NEUMAIER ''AS IS'' AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL  ARNOLD NEUMAIER BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobnan.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = snobnan(fnan,f,near,inew)
% replaces the function values NaN of a set of points by a value 
% determined by their nearest neighbors with finite function values, 
% with a safeguard for the case that all neighbors have function value
% NaN
%
% Input:
% fnan	  	vector containing the pointers to the points where the
%               function value could not be obtained
% f	  	f(:,1) set of available function values
%		f(:,2) their uncertainty/variation
% near(j,:)	vector pointing to the nearest neighbors of point j
% inew          vector pointing to the new boxes and boxes whose nearest
%               neighbors have changed
%
% Output:
% f		updated version of f
%
function f = snobnan(fnan,f,near,inew)
lnn = size(near,2);
notnan = 1:size(f,1);
notnan(fnan) = [];
[fmx,imax] = max(f(notnan,1));
fmn = min(f(notnan,1));
dfmax = f(imax,2);
for j=1:length(fnan)
  l = fnan(j);
  if ~isempty(find(inew==l))
% a substitute function value is only computed for new points and for 
% points whose function values have changed
    ind = near(l,:); 
% eliminate neighbors with function value NaN
    ind1 = []; 
    for i=1:length(ind)
      if ~isempty(find(fnan==ind(i))),ind1 = [ind1 i]; end
    end
    if ~isempty(ind1)
      ind(ind1) = [];
    end
    if isempty(ind)
      f(l,1) = fmx + 1.e-3*(fmx-fmn);
      f(l,2) = dfmax;
    else
      [fmax1,k] = max(f(ind,1));
      fmin1 = min(f(ind,1));
      f(l,1) = fmax1 + 1.e-3*(fmax1-fmin1);
      f(l,2) = f(ind(k),2);
    end
  end
end

