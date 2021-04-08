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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rsort.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,w,cdfx,dof]=rsort(x,w);
% sort x in increasing order
% remove multiple entries, and adapt weights w accordingly
% x and w must both be a row or a column
% default input weights are w=1
% 
% if nargout>2, the weighted empirical cdf is computed at x
% dof = length(x) at input is the original number of degrees of freedom
%
function [x,w,cdfx,dof]=rsort(x,w);

if nargin==1, w=ones(size(x)); end;

debug=0;
if debug>0, 
  save rsort.debug x w
elseif debug<0,
  load rsort.debug; 
  disp('debug mode');
end;


% turn into rows and sort
col=(size(x,2)==1);
if col, x=x';w=w'; end;
[x,ind]=sort(x);
w=w(ind);

% remove repetitions
n=length(x);
ind=find([x(2:n),inf]~=x);
nn=length(ind);
x=x(ind);
w(1)=sum(w(1:ind(1)));
for i=2:nn,
  w(i)=sum(w(ind(i-1)+1:ind(i)));
end;
w(nn+1:n)=[];

% restore original shape
if col, x=x';w=w'; end;

if nargout<3, return; end; 

% get cumulative sum of weights
cdfx=w;
for i=2:nn,
  cdfx(i)=cdfx(i-1)+w(i);
end;

% adjust for jumps and normalize
cdfx=(cdfx-0.5*w)/cdfx(nn);
dof=n;



return 

%%%%%%%% test documentation %%%%%%%%%
% test rsort.m (outside of rsort.m)

n=200;

u=[1 1 1 5 9 13 27 66 78 85 88 93]';

figure(1);clf
while 1,
  x=u(fix(10*rand(n,1))+1);
  [x,w,cdfx,dof]=rsort(x);
  [x,w]
  plot(x,cdfx)
  input('next>');
end;
