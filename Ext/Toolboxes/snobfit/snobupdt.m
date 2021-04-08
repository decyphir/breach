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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobupdt.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xl,xu,x,f,nsplit,small,near,d,np,t,inew,fnan,u,v] = 
% snobupdt(xl,xu,x,f,nsplit,small,near,d,np,t,xnew,fnew,fnan,u,v,u1,v1,
% dx)
% updates the box parameters when a set of new points and their function
% values are added, i.e., the boxes containing more than one point are
% split and the nearest neighbors are computed or updated
%
% Input:
% xl,xu		rows contain lower and upper bounds of the old boxes
% x		rows contain the old points
% f		f(j,:) contains the function value at x(j,:), its
%		variation and other parameters
% nsplit	nsplit(j,i) number of times box j has been split along
%		the ith coordinate
% small		small(j) is an integer-valued logarithmic volume measure
%		of box j
% near		near(j,:) is a vector pointing to the nearest 
%		neighbors of x(j,:)
% d		d(j) is the maximal distance between x(j,:) and one of
%		its neighbors
% np	        np(j) is the number of times the function value of 
%		x(j,:) has been measured
% t		t(j) is np(j) times the variance of the function values
%		measured for the point x(j,:)
% xnew		rows contain new points
% fnew		new function values and their variations
%		fnew(j,1) = f(xnew(j,:)), fnew(j,2) = df(xnew(j,:))
% fnan          pointer to all old points where the function value could
%               not be obtained
% u,v		box bounds
% u1,v1         box in which the new points are to be generated
% dx            resolution vector
%
% Output:
% xl,xu		updated version of xl,xu (including new boxes)
% x		updated version of x
% f		updated version of f
% nsplit	updated version of nsplit
% small		updated version of small
% near		updated version of near
% d		updated version of d
% np		updated version of np
% t		updated version of t
% inew		pointer pointing to the new boxes and boxes whose 
%               nearest neighbors have changed
% fnan          possibly updated version of fnan (if a function value
%               was found for a point in the new iteration)
% u,v		possibly updated box bounds such that all new points
% 		are in the box
%
function [xl,xu,x,f,nsplit,small,near,d,np,t,inew,fnan,u,v] = snobupdt(xl,xu,x,f,nsplit,small,near,d,np,t,xnew,fnew,fnan,u,v,u1,v1,dx)
n = length(u);  %dimension of the problem
nneigh = n+5;
nxold = size(x,1);	% number of points from the previous iteration
nxnew = size(xnew,1);
inew = [];
if nxold
% if any of the new points are already among the old points, they are
% thrown away and the function value and its uncertainty are updated
  del = [];
  for j=1:nxnew
    i = find(sum(abs(ones(nxold,1)*xnew(j,:)-x),2)==0);
    if ~isempty(i)
      if isempty(find(fnan==i)) & isfinite(f(i,1)) % point i had finite
                                                   % function value
        if ~isnan(fnew(j,1))
          np(i) = np(i) + 1;
          delta = fnew(j,1) - f(i,1);
          f(i,1) = f(i,1) + delta/np(i);
          t(i) = t(i) + delta*(fnew(j,1)-f(i,1));
          f(i,2) = sqrt(f(i,2)^2+(delta*(fnew(j,1)-f(i,1))+fnew(j,2)^2-f(i,2)^2)/np(i));
          inew = [inew i];
        end
        del = [del j];
      else % point i had NaN function value
        if ~isnan(fnew(j,1))
          f(i,1:2) = fnew(j,1:2);
          inew = [inew i];
          ii = find(fnan==i);
          if ~isempty(ii), fnan(ii) = []; end
        end
        del = [del j];
      end
    end
  end
  xnew(del,:) = [];
  fnew(del,:) = [];
  nxnew = size(xnew,1);
  if ~nxnew, inew = sort(inew); return, end
end
[xnew,fnew,npnew,tnew] = snobinput(xnew,fnew);
nxnew = size(xnew,1);
nx = nxold + nxnew;	% current number of points
if sum(min([xnew; u],[],1)<u) | sum(max([xnew;v],[],1)>v) ...
   | sum(min(u,u1)<u) | sum(max(v,v1)>v)
  [xl,xu,small,u,v] = snobnewb(xnew,xl,xu,small,u,v,u1,v1);
end
x = [x; xnew];
inew = [inew nxold+1:nx];
f(nxold+1:nx,1:2) = fnew;
np = [np npnew];
t = [t tnew];
if ~nxold
  [xl,xu,x,f,nsplit,small] = snobsplit(x,f,u,v,zeros(1,n),u,v);
else
  for j=1:nxnew
    xx = ones(nxold,1)*xnew(j,:);
    ind = find(sum(xl<=xx&xx<=xu,2)==n);
    [minsmall,ismall] = min(small(ind));
    par(j) = ind(ismall);
  end
  par1 = rsort(par);
  inew = [inew par1];
  for l=1:length(par1)
    j = par1(l);
    ind = find(par==j);
    ind = ind + nxold;
    [xl0,xu0,x0,f0,nsplit0,small0] = snobsplit...
    (x([j ind],:),f([j ind],1:2),xl(j,:),xu(j,:),nsplit(j,:),u,v);
    nxj = length(ind)+1;     % number of points in box [xl(j,:),xu(j,:)]
    k = find(sum(x0==ones(nxj,1)*x(j,:),2)==n);
    xl(j,:) = xl0(k,:);
    xu(j,:) = xu0(k,:);
    nsplit(j,:) = nsplit0(k,:);
    small(j) = small0(k);
    for k=1:nxj-1
      k1 = ind(k);
      k2 = find(sum(x0==ones(nxj,1)*x(k1,:),2)==n);
      if length(k2) > 1
        [msmall,k3] = min(small(k2));
        k2 = k2(k3);
      end
      xl(k1,:) = xl0(k2,:);
      xu(k1,:) = xu0(k2,:);
      nsplit(k1,:) = nsplit0(k2,:);
      small(k1) = small0(k2);
    end
  end
end
notnan = 1:nx; 
fnan1 = find(isnan(f(:,1)));
notnan([fnan fnan1]) = [];
if ~isempty(notnan)
  fmn = min(f(notnan,1));
  fmx = max(f(notnan,1));
else
  fmn = 1;
  fmx = 0;
end
if nx >= nneigh+1 & fmn < fmx
  for j=nxold+1:nx
    [near(j,:),d(j)] = snobnn(x(j,:),x,nneigh,dx);
  end
  for j=1:nxold
    if min(sqrt(sum((ones(nxnew,1)*x(j,:)-xnew).^2,2))) < d(j)
      [near(j,:),d(j)] = snobnn(x(j,:),x,nneigh,dx);
      inew = [inew j];
    end
  end
  inew = sort(inew);
else
  near = [];
  d = Inf*ones(1,nx);
end
