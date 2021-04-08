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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobnn.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [near,d] = snobnn(x0,x,m,dx)
% computes a safeguarded set of m nearest neighbors to a point x0
% for each coordinate, the nearest point differing from x0 in the ith
% coordinate by at least 1 is chosen
% 
% Input:
% x0	point for which the neighbors are to be found
% x	the rows are the points from which the neighbors are to be
%	chosen (possibly x0 is among these points)
% m	number of neighbors to be found
% dx    resolution vector
%
% Output:
% near	vector pointing to the nearest neighbors (i.e. to the rows
%	of x)
% d	maximal distance between x0 and a point of the set of m
% 	nearest neighbors
% 
function [near,d] = snobnn(x0,x,m,dx)
n = length(x0);  % dimension of the problem
d = sqrt(sum((x-ones(size(x,1),1)*x0).^2,2));
[d1,ind] = sort(d);
if ~d1(1), ind(1) = []; end	% eliminate x0 if it is in the set
near = [];
for i=1:n
  j = min(find(abs(x(ind,i)-x0(i))-dx(i)>=0));
  near = [near ind(j)];
  ind(j) = [];
end
j = 1;
while length(near) < m
  if isempty(find(near==ind(j))) & max(abs(x0-x(ind(j),:))-dx)>= 0
    near = [near ind(j)]; 
  end
  j = j + 1;
end
[d,ind]=sort(d(near));
near = near(ind);
d = max(d);
