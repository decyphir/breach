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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobqfit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [y,f1] = snobqfit(j,x,f,near,dx,u,v)
% a quadratic model around the best point is fitted and minimized with
% minq over a trust region
%
% Input:
% j     index of the best point
% x     the rows contain the points where the function has been
%       evaluated
% f     corresponding function values, i.e., f(i) = f(x(i,:))
% near  near(i,:) is a vector containing the indices of the nearest
%       neighbors of the point x(i,:)
% dx    resolution vector, i.e. the ith coordinate of a point to be
%       generated is an integer-valued multiple of dx(i) 
% u,v   the points are to be generated in [u,v]
%
% Output:
% y     minimizer of the quadratic model around the best point
% f1    its estimated function value
%
function [y,f1] = snobqfit1(j,x,f,near,dx,u,v)
n = size(x,2);	% dimension of the problem
K = min(size(x,1)-1,n*(n+3));
nneigh = size(near,2);
x0 = x(j,:);
f0 = f(j);
distance = sum((x-ones(size(x,1),1)*x0).^2,2);
[dd,ind] = sort(distance);
ind = ind(2:K+1);
d = max(abs(x(near(j,:),:)-ones(nneigh,1)*x0),[],1);
d = max(d,dx);
S = x(ind,:) - ones(K,1)*x0;
R = triu(qr(S,0));
R = R(1:n,:);
L = inv(R)';
sc=sum((S*L').^2,2).^(3/2);
b = (f(ind)-f0)./sc;
A = [x(ind,:) - ones(K,1)*x0 0.5*(x(ind,:) - ones(K,1)*x0).^2];
for i=1:n-1
  B = (x(ind,i)-x0(i))*ones(1,n-i);
  A = [A B.*(x(ind,i+1:n)-ones(K,1)*x0(i+1:n))];
end
A = A./(sc*ones(1,size(A,2)));
y = A\b;
for i=1:n
  G(i,i) = y(n+i);
end
l = 2*n+1;
for i = 1:n-1
  for j=i+1:n
    G(i,j) = y(l);
    G(j,i) = y(l);
    l = l + 1;
  end
end
g = y(1:n);
[y,f1] = minq(f0-x0*g+0.5*x0*G*x0',g-G*x0',G,max(x0'-d',u'),min(x0'+d',v'),0);
y = snobround(y',u,v,dx);
nc = 0;
while min(max(abs(x-ones(size(x,1),1)*y)-ones(size(x,1),1)*dx,[],2))<-eps & nc<10
  u1 = max(x0-d,u);
  v1 = min(x0+d,v);
  y = u1+rand(1,n).*(v1-u1);
  y = snobround(y,u,v,dx);%disp(y)
  nc = nc + 1;
end
f1 = f0+(y-x0)*g+0.5*(y-x0)*G*(y-x0)';
