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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snob5.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x1 = snob5(x,u,v,dx,nreq)
% generates nreq points of class 5 in [u,v]
%
% Input:
% x       the rows are the points already chosen by Snobfit (can be 
%         empty)
% u,v     bounds of the box in which the points are to be generated
% dx      resolution vector, i.e. the ith coordinate of a point to be
%         generated is an integer-valued multiple of dx(i) 
% nreq    number of points to be generated
%
% Output:
% x1      the rows are the requested points
%         x1 is of dimension nreq x n (n = dimension of the problem)
%
function x1 = snob5(x,u,v,dx,nreq)
n = length(u);   % dimension of the problem
nx = size(x,1);
nx1 = 100*nreq;
xnew = ones(nx1,1)*u+rand(nx1,n).*(ones(nx1,1)*(v-u));
if nx
  for j=1:nx1
    xnew(j,:) = snobround(xnew(j,:),u,v,dx);
    d(j) = min(sum((x-ones(nx,1)*xnew(j,:)).^2,2));
  end
  ind = find(~d);
  xnew(ind,:) = [];
  d(ind) = [];
  x1 = [];
  nx1 = size(xnew,1);
  if size(d,2) > 1, d = d'; end
else
  x1 = xnew(1,:);
  xnew(1,:) = [];
  nx1 = nx1-1;
  d = sum((xnew-ones(nx1,1)*x1).^2,2);
  nreq = nreq -1;
end
for j = 1:nreq
  [dmax,i] = max(d);
  y = xnew(i,:);
  x1 = [x1; y];
  xnew(i,:) = [];
  d(i) = [];
  nx1 = nx1 - 1;
  d1 = sum((xnew-ones(nx1,1)*y).^2,2);
  d = min(d,d1);
end
