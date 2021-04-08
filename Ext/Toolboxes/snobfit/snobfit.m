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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobfit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [request,xbest,fbest] = snobfit(file,x,f,params,dx)
% minimization of a function over a box in R^n
%
% Input:
% file		name of file for input and output
%		if nargin < 5, the program continues a previous run and
%		reads from file.mat
%		the output is (again) stored in file.mat
% x		the rows are a set of new points entering the 
%		optimization algorithm together with their function 
%		values  
% f		matrix containing the corresponding function values
%		and their uncertainties, i.e., f(j,1) = f(x(j,:)) 
% 		and f(j,2) = df(x(j,:))
%               a value f(j,2)<=0 indicates that the corresponding
%               uncertainty is not known, and the program resets it to
%               sqrt(eps)
% params	structure variable defining the box [u,v] in which the
%		points are to be generated, the number nreq of 
%		points to be generated and the probability p that a 
%               point of type 4 is generated
%		params = struct('bounds',{u,v},'nreq',nreq,'p',p)
% dx		only used for the definition of a new problem (when
%		the program should continue from the values stored in
%		file.mat, the call should have only 4 input parameters!)
%	        n-vector (n = dimension of the problem) of minimal
%               steps, i.e., two points are considered to be different 
%               if they differ by at least dx(i) in at least one 
%               coordinate i
%
% Output:
% request	nreq x (n+3)-matrix
%		request(j,1:n) is the jth newly generated point, 
%		request(j,n+1) is its estimated function value and
%		request(j,n+3) indicates for which reason the point
%		request(j,1:n) has been generated 
%		request(j,n+3) = 1 best prediction
%		               = 2 putative local minimizer
%		               = 3 alternative good point
%		               = 4 explore empty region
%                              = 5 to fill up the required number of
%                              function values if too little points of
%                              the other classes are found
% xbest		current best point
% fbest		current best function value (i.e. function value at 
%		xbest)
%
% Uses the following m-files (directly or indirectly):
% minq.m and its subprograms
% rsort.m
% snob5.m
% snobinput.m
% snoblocf.m
% snoblp.m
% snobnan.m
% snobnewb.m
% snobnn.m
% snobpoint.m
% snobqfit.m
% snobqmin.m
% snobround.m
% snobsplit.m
% snobupdt.m
% snobwarn.m
%
function [request,xbest,fbest] = snobfit(file,x,f,params,dx)

if isempty(f), f=zeros(0,2); end;
ind = find(f(:,2)<=0);
if ~isempty(ind)
  f(ind,2) = sqrt(eps);
end
rho = 0.5*(sqrt(5)-1);	% golden section number
u1 = params(1).bounds;
v1 = params(2).bounds;
nreq = params(1).nreq;
p = params(1).p;
n = length(u1);	% dimension of the problem
nneigh = n+5;   % number of nearest neighbors
if size(u1,1) > 1, u1 = u1'; end
if size(v1,1) > 1, v1 = v1'; end
dy = 0.1*(v1-u1);% defines the vector of minimal distances between two 
                 % points suggested in a single call to Snobfit
if nargin > 4	 % a new job is started
  if any(dx<=0), error('dx should only contain positive entries'), end
  if size(dx,1) > 1, dx = dx'; end
  if ~isempty(x)
    u = min([x(:,1:n);u1]);
    v = max([x(:,1:n);v1]);
  else
    u = u1;
    v = v1;
  end                     
  [x,f,np,t] = snobinput(x,f); % throw out duplicates among the points
                               % and compute mean function value and
                               % deviation 
  if size(x,1) 
    [xl,xu,x,f,nsplit,small] = snobsplit(x,f,u,v,zeros(1,n),u,v);
    d = Inf*ones(1,size(x,1));
  else
    xl = [];
    xu = [];
    nsplit = [];
    small = [];
  end
  notnan = find(isfinite(f(:,1)));
  if ~isempty(notnan)
    fmn = min(f(notnan,1));
    fmx = max(f(notnan,1));
  else
    fmn = 1;
    fmx = 0;
  end
  if size(x,1) >= nneigh+1 & fmn < fmx
    inew = 1:size(x,1);
    for j=inew
      [near(j,:),d(j)] = snobnn(x(j,:),x,nneigh,dx);
    end
    fnan = find(isnan(f(:,1)));
    if ~isempty(fnan)
      f = snobnan(fnan,f,near,inew);
    end
    for j=inew
       [y(j,:),f(j,3),g(j,:),sigma(j)] = snoblocf(j,x,f,near,dx,u,v);
    end
    [fbest,jbest] = min(f(:,1));
    xbest = x(jbest,:);
  else
    fnan = [];
    near = [];
    d = Inf*ones(1,size(x,1));
    x1 = snob5(x,u1,v1,dx,nreq);
    request = [x1 NaN*ones(nreq,1) 5*ones(nreq,1)];
    if ~isempty(x)
      [fbest,jbest] = min(f(:,1));
      xbest = x(jbest,:);
    else
      xbest = NaN*ones(1,n);
      fbest = Inf;
    end
    if size(request,1) < nreq, snobwarn, end
    clear xnew fnew u1 v1 nreq
    save(file)
    return
  end 
else
  xnew = x;
  fnew = f;
  load(file)
  nx = size(xnew,1);
  oldxbest = xbest;
  [xl,xu,x,f,nsplit,small,near,d,np,t,inew,fnan,u,v] = snobupdt(xl,xu,x,f,nsplit,small,near,d,np,t,xnew,fnew,fnan,u,v,u1,v1,dx);
  if ~isempty(near) 
    ind = find(isnan(f(:,1)));
    fnan = [fnan; ind];
    if ~isempty(fnan)
      f = snobnan(fnan,f,near,inew);
    end
    [fbest,jbest] = min(f(:,1));
    xbest = x(jbest,:);
    for j=inew
      [y(j,:),f(j,3),g(j,:),sigma(j)] = snoblocf(j,x,f(:,1:2),near,dx,u,v);
    end
  else
    x1 = snob5(x,u1,v1,dx,nreq);
    request = [x1 NaN*ones(nreq,1) 5*ones(nreq,1)];
    if ~isempty(x)
      [fbest,ibest] = min(f(:,1));
      xbest = x(ibest,:);
    else
      xbest = [];
      fbest = Inf;
    end
    if size(request,1) < nreq, snobwarn, end
    clear xnew fnew u1 v1 nreq
    save(file)
    return
  end
end

sx = size(x,1);
request = [];
ind = find(sum(xl<=ones(sx,1)*v1&xu>=ones(sx,1)*u1,2)==n);
[minsmall,k] = min(small(ind));
maxsmall=max(small(ind));
m1 = floor((maxsmall-minsmall)/3);
k = find(small(ind)==minsmall);
k = ind(k);
[fsort,j] = sort(f(k,1));
k = k(j);
isplit = k(1);

if sum(u1<=xbest&xbest<=v1) == n
  [z,f1] = snobqfit(jbest,x,f(:,1),near,dx,u1,v1);
  z = snobround(z,u1,v1,dx);
  zz = ones(sx,1)*z;
  j = find(sum(xl<=zz&zz<=xu,2)==n);
  if length(j) > 1
    [msmall,j1] = min(small(j));
    j = j(j1);
  end
  if min(max(abs(x-ones(sx,1)*z)-ones(sx,1)*dx,[],2)) >= -eps
    dmax = max((xu(j,:)-xl(j,:))./(v-u));
    dmin = min((xu(j,:)-xl(j,:))./(v-u));
    if dmin <= 0.05*dmax
      isplit = [isplit j];
    else 
      request = [request; z f1 1];
    end
  end
else
  [fbes,jbes] = min(f(ind,1));
  jbes = ind(jbes);
  xbes = x(jbes,:);
  [z,f1] = snobqfit(jbes,x,f(:,1),near,dx,u1,v1);
  z = snobround(z,u1,v1,dx);
  zz = ones(size(x,1),1)*z;
  j = find(sum(xl<=zz&zz<=xu,2)==n);
  if length(j) > 1
    [msmall,j1] = min(small(j));
    j = j(j1);
  end
  if min(max(abs(x-ones(sx,1)*z)-ones(sx,1)*dx,[],2)) >= -eps
    dmax = max((xu(j,:)-xl(j,:))./(v-u));
    dmin = min((xu(j,:)-xl(j,:))./(v-u));
    if dmin <= 0.05*dmax
      isplit = [isplit j];
    else
      request = [request; z f1 1];
    end
  end
end

if size(request,1) < nreq 
  globloc = nreq-size(request,1);
  glob1 = globloc*p;
  glob2 = floor(glob1);
  if rand < glob1 - glob2;
    glob = glob2 + 1;
  else
    glob = glob2;
  end
  loc = globloc - glob;
  if loc
    [local,nlocal] = snoblp(f(:,1),near,ind);
    [fsort,k] = sort(f(local,3));
    j = 1;
    sreq = size(request,1);
    while sreq < nreq-glob & j <= length(local)
      l0 = local(k(j));
      y1 = snobround(y(l0,:),u1,v1,dx);
      yy = ones(size(x,1),1)*y1;
      l = find(sum(xl<=yy&yy<=xu,2)==n);
      if length(l) > 1
        [msmall,j1] = min(small(l));
        l = l(j1);
      end
      dmax = max((xu(l,:)-xl(l,:))./(v-u));
      dmin = min((xu(l,:)-xl(l,:))./(v-u));
      if dmin <= 0.05*dmax
        isplit = [isplit l];
        j = j + 1;
        continue
      end
      if max(abs(y1-x(l,:))-dx)>=-eps & (~sreq | min(max(abs(request(:,1:n)-ones(sreq,1)*y1)-ones(sreq,1)*max(dy,dx),[],2))>=-eps) 
        if sum(y1==y(l0,:)) < n
          D = f(l0,2)./dx.^2;
          f1 = f(l0,1)+g(l0,:)*(y1-x(l0,:))'+sigma(l0)*((y1-x(l0,:))*diag(D)*(y1-x(l0,:))'+f(l0,2));
        else
          f1 = f(l0,3);
        end
        request = [request; y1 f1 2]; 
      end
      sreq = size(request,1);
      j = j + 1;      
    end
    if sreq < nreq-glob
      [fsort,k] = sort(f(nlocal,3));
    end
    j = 1;
    while sreq < nreq-glob & j<= length(nlocal) 
      l0 = nlocal(k(j));
      y1 = snobround(y(l0,:),u1,v1,dx);
      yy = ones(size(x,1),1)*y1;
      l = find(sum(xl<=yy&yy<=xu,2)==n);
      if length(l) > 1
        [msmall,j1] = min(small(l));
        l = l(j1);
      end
      dmax = max((xu(l,:)-xl(l,:))./(v-u));
      dmin = min((xu(l,:)-xl(l,:))./(v-u));
      if dmin <= 0.05*dmax
        isplit = [isplit l];
        j = j + 1;
        continue
      end
      if max(abs(y1-x(l,:))-dx)>=-eps & (~sreq | min(max(abs(request(:,1:n)-ones(sreq,1)*y1)-ones(sreq,1)*max(dy,dx),[],2))>=-eps) 
        if sum(y1==y(l0,:)) < n
          D = f(l0,2)./dx.^2;
          f1 = f(l0,1)+g(l0,:)*(y1-x(l0,:))'+sigma(l0)*((y1-x(l0,:))*diag(D)*(y1-x(l0,:))'+f(l0,2));
        else
          f1 = f(l0,3);
        end 
        request = [request; y1 f1 3];
      end
      sreq = size(request,1);
      j = j + 1;      
    end
  end
end
sreq = size(request,1);
for l=isplit
  jj = find(ind==l);
  ind(jj) = [];
  [y1,f1] = snobpoint(x(l,:),xl(l,:),xu(l,:),f(l,1:2),g(l,:),sigma(l),u1,v1,dx);
  if max(abs(y1-x(l,:))-dx)>=-eps & (~sreq | min(max(abs(request(:,1:n)-ones(sreq,1)*y1)-ones(sreq,1)*dx,[],2))>=-eps) 
    request = [request; y1 f1 4];
  end
  sreq = size(request,1);
  if sreq == nreq, break, end
end
first = 1;
while sreq < nreq & ~isempty(ind) %& ~isempty(find(small(ind)<=minsmall+m1))
  for m=0:m1
    if first == 1 
      first = 0;
      continue 
    end
    m = 0;
    k = find(small(ind)==minsmall+m);
    while isempty(k)
      m = m+1;
      k = find(small(ind)==minsmall+m);
    end
    if ~isempty(k)
      k = ind(k);
      [fsort,j] = sort(f(k,1));
      k = k(j);
      l = k(1);
      jj=find(ind==l);
      ind(jj)=[];
      [y1,f1] = snobpoint(x(l,:),xl(l,:),xu(l,:),f(l,1:2),g(l,:),sigma(l),u1,v1,dx);
      if  max(abs(y1-x(l,:))-dx)>=-eps & (~sreq | min(max(abs(request(:,1:n)-ones(sreq,1)*y1)-ones(sreq,1)*max(dy,dx),[],2))>=-eps)  
        request = [request; y1 f1 4];
      end
      sreq = size(request,1);
      if sreq == nreq, break,end
    end
    m = 0;
  end
end
if size(request,1) < nreq
  x1 = snob5([x; request(:,1:n)],u1,v1,dx,nreq-size(request,1));
  nx = size(x,1);
  for j=1:size(x1,1)
    i = find(sum(xl<=ones(nx,1)*x1(j,:)&ones(nx,1)*x1(j,:)<=xu,2)==n);
    if length(i) > 1
      [minv,i1] = min(small(i));
      i = i(i1);
    end
    D = f(i,2)./dx.^2;
    f1 = f(i,1) + (x1(j,:)-x(i,:))*g(i,:)' + sigma(i)*((x1(j,:)-x(i,:))*diag(D)*(x1(j,:)-x(i,:))'+f(i,2));
    request = [request; x1(j,:) f1 5];
  end
end
if size(request,1) < nreq, snobwarn, end
clear xnew fnew u1 v1 nreq p
save(file)


