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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobsplit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xl,xu,x,f,nsplit,small] = snobsplit(x,f,xl0,xu0,nspl,u,v)
% splits a box [xl0,xu0] contained in a bigger box [u,v] such that
% each of the resulting boxes contains just one point of a given set of
% points
%
% Input:
% x	the rows are a set of points
% f	f(j,:) contains the function value, its variation and possibly
%	other parameters belonging to point x(j,:)
% xl0	vector of lower bounds of the box
% xu0	vector of upper bounds of the box
% nspl	nspl(i) is the number of splits the box [xl0,xu0] has already
%	undergone along the ith coordinate
% u,v	bounds of the original problem	
%	[xl0,xu0] is contained in [u,v]
%
% Output:
% xl	xl(j,:) is the lower bound of box j
% xu	xi(j,:) is the upper bound of box j
% x      	x(j,:) is the point contained in box j
% f	f(j,:) contains the function value at x(j,:), its uncertainty
%	etc. 
% nsplit nsplit(j,i) is the number of times box j has been split in the
%	ith coordinate
% small	small(j) = integer-valued logarithmic volume measure of box j
%
function [xl,xu,x,f,nsplit,small] = snobsplit(x,f,xl0,xu0,nspl,u,v)

n = size(x,2);
if size(xl0,1) > 1
  xl0 = xl0';
  xu0 = xu0';
end
if nargin < 5
  nspl = zeros(1,n);
end
if nargin < 6
  u = xl0;
  v = xu0;
end
if size(x,1) == 1
  xl = xl0;
  xu = xu0;
  nsplit = nspl;
  small = -sum(round(log2((xu-xl)./(v-u))));
  return
elseif size(x,1) == 2
  [dmax,i] = max(abs(x(1,:)-x(2,:))./(v-u));
  ymid = 0.5*(x(1,i)+x(2,i));
  xl(1,:) = xl0;
  xu(1,:) = xu0;
  xl(2,:) = xl0;
  xu(2,:) = xu0;
  if x(1,i) < x(2,i)
    xu(1,i) = ymid;
    xl(2,i) = ymid;
  else
    xl(1,i) = ymid;
    xu(2,i) = ymid;
  end
  nsplit(1,:) = nspl;
  nsplit(1,i) = nsplit(1,i) + 1;
  nsplit(2,:) = nsplit(1,:);
  small = -sum(round(log2((xu-xl)./(ones(2,1)*(v-u)))),2)';
  return
end
for i=1:n
  var(i) = std(x(:,i)/(v(i)-u(i)));
end 
xx = sort(x,1);
dd = xx(2:size(xx,1),:)-xx(1:size(xx,1)-1,:);
[mvar,i] = max(var);
y = rsort(x(:,i));
d = y(2:length(y)) - y(1:length(y)-1);
ld = length(d);
ii=1:ld;
[dmax,j] = max(d(ii));
j = ii(j);
ymid = 0.5*(y(j)+y(j+1));
ind1 = find(x(:,i)<ymid);
ind2 = find(x(:,i)>ymid);
xl(1,:) = xl0;
xu(1,:) = xu0;
xu(1,i) = ymid;
xl(2,:) = xl0;
xl(2,i) = ymid;
xu(2,:) = xu0;
nsplit(1,:) = nspl;
nsplit(1,i) = nsplit(1,i) + 1;
nsplit(2,:) = nsplit(1,:);
npoint(1) = length(ind1);
npoint(2) = length(ind2);
ind(1,1:length(ind1)) = ind1';
ind(2,1:length(ind2)) = ind2';
nboxes = 2;
[maxpoint,j] = max(npoint);
while  maxpoint > 1
  ind0 = ind(j,find(ind(j,:)));
  for i=1:n
    var(i) = std(x(ind0,i)/(v(i)-u(i)));
  end
  [maxvar,i] = max(var);
  y = rsort(x(ind0,i));
  d = y(2:length(y)) - y(1:length(y)-1);
  ld = length(d);
  ii=1:ld;
  [dmax,k] = max(d(ii));
  k = ii(k);
  ymid = 0.5*(y(k)+y(k+1));
  ind1 = find(x(ind0,i)<ymid);
  ind2 = find(x(ind0,i)>ymid);
  ind1 = ind0(ind1);
  ind2 = ind0(ind2);
  nboxes = nboxes + 1;
  xl(nboxes,:) = xl(j,:);
  xu(nboxes,:) = xu(j,:);
  xu(j,i) = ymid;
  xl(nboxes,i) = ymid;
  nsplit(j,i) = nsplit(j,i) + 1;
  nsplit(nboxes,:) = nsplit(j,:);
  npoint(j) = length(ind1);
  npoint(nboxes) = length(ind2);
  ind(j,1:size(ind,2)) = 0;
  ind(j,1:length(ind1)) = ind1;
  ind(nboxes,1:length(ind2)) = ind2;
  if ~sum(ind(:,size(ind,2))), ind(:,size(ind,2)) = []; end
  [maxpoint,j] = max(npoint);
end
x = x(ind(:,1),:);
f = f(ind(:,1),:);
small = -sum(round(log2((xu-xl)./(ones(size(x,1),1)*(v-u)))),2)';

