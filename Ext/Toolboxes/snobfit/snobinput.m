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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobinput.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,np,t] = snobinput(x,f)
% checks whether there are any duplicates among the points given by the
% rows of x, throws away the duplicates and computes their average
% function values and an estimated uncertainty
%
% Input:
% x	the rows of x are a set of points
% f	f(j,1) is the function value of x(j,:) and f(j,2) is its
%	uncertainty
%
% Output:
% x	updated version of x (possibly some points have been deleted)
% f	updated version of f (f(j,1) is the average function value and
%	f(j,2) the estimated uncertainty pertaining to x(j,:))  
% np	np(j) is the number of times the row x(j,:) appeared in the
%	input version of x
% t	t(j) is np(j) times the variance of the function values measured
%	for point x(j,:)
%
function [x,f,np,t] = snobinput(x,f)
sx = size(x,1);
n = size(x,2);
i = 1;
if isempty(x)
  np = [];
  t = [];
end
while i <= sx
  j = i + 1;
  ind = [];
  while j <= sx
    if sum(x(i,:)==x(j,:)) == n
      ind = [ind j];
    end
    j = j + 1;
  end
  if ~isempty(ind)
    ind = [i ind];
    ind1 = find(isnan(f(ind,1)));
    if length(ind1) < length(ind)
      ind(ind1) = [];
      np(i) = length(ind);
      fbar = sum(f(ind,1))/np(i);
      t(i) = sum((f(ind,1)-fbar).^2);
      f(i,1) = fbar;
      f(i,2) = sqrt((sum(f(i,2).^2)+t(i))/np(i));
    else
      np(i) = 1;
      t(i) = 0;
    end
    x(ind(2:end),:) = [];
    f(ind(2:end),:) = [];
    sx = size(x,1);
  else
    np(i) = 1;
    t(i) = 0;
  end 
  i = i + 1;
end






