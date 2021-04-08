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
%%%%%%%%%%%%%%%%%%%%%%%%%% snobsofttest.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file (m-script) for applying the soft optimality theorem and SNOBFIT
% to some Hock-Schittkowski problems
clear; clear global; clear mex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beginning of data to be adapted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nprob = 18;     % problem number (nprob in {18,19,23,30,31,34,36,37,41,
                % 53,60,65,66,71,73,74})
fun = ['hsf' int2str(nprob)];  % name of objective function
Fun = ['hsF' int2str(nprob)];  % name of constraint function
file = 'hs';    % the intermediate data (after each call to SNOBFIT)
                % are stored in hs.mat
[u,v,F1,F2,x0,fglob,xglob] = hsdata(nprob); 
% bounds for the variables and the constraints, default starting point,
% optimal function value and global minimizer
n = length(u);  % dimension of the problem
m = length(F1); % number of constraints 
ncall = 10000;   % limit on the number of function calls
npoint = n+6;   % the standard starting point and npoint-1 randomly
                % generated points are used as input for the initial
                % call to SNOBFIT
nreq = n+6;     % number of points to be generated in the continuation
                % calls to SNOBFIT
dx = (v-u)'*1.e-5; % resolution vector
p = 0.5;        % fraction of points of class 4 among the points of 
                % classes 2 to 4
sigma = 0.05*[25; 25]; % choice of sigma in the merit function
% (column vector)
prt = 0;        % print level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of data to be adapted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = struct('bounds',{u,v},'nreq',nreq,'p',p); % input structure
x1 = [];
f = [];
F = [];
change = 0;
% generation of npoint-1 random points in [u,v]
x = rand(npoint-1,n);
x = x*diag(v-u)+ones(npoint-1,1)*u';
x = [x0';x];
% computation of the objective function values and constraints for the 
% starting points
f0 = Inf; % initialization of f0
for j=1:npoint
  x(j,:) = snobround(x(j,:),u',v',dx);
  ff = feval(fun,x(j,:));
  FF = feval(Fun,x(j,:));
  if size(FF,1) > 1
    FF = FF';
  end
  f = [f; ff];
  F = [F; FF];
  if (sum(F1'<=FF&FF<=F2')==n)
    f0 = min(f0,ff); % f0 is the smallest function value among the
                     % feasible ones of the initial points
  end 
end
x1 = [x1; x];
if isinf(f0) 
  f0 = 2*max(f)-min(f);
end
Delta = median(abs(f-f0));
% computation of the merit function for the starting points
for j=1:npoint
  fm(j,1) = softmerit(f(j),F(j,:),F1,F2,f0,Delta,sigma);
end
fm(:,2) = sqrt(eps)*ones(npoint,1); % set uncertainty of the merit 
                                    % function
ncall0 = npoint; % function call counter
while ncall0 < ncall  % repeat until the limit on function calls is 
                      % reached
  if ncall0 == npoint % initial call
    [request,xbest,fbest] = snobfit(file,x,fm,params,dx);
    ncall0,xbest,fbest
  else                % continuation call
    [request,xbest,fbest] = snobfit(file,x,fm,params);
  end
  if prt, request, end         % shows the suggested points, estimated
                               % merit function values and class (rowwise)
  clear x; clear fm
  % compute merit function at the suggested points
  for j=1:size(request,1)
    x(j,:) = request(j,1:n);
    x1 = [x1; x(j,:)];
    ff = feval(fun,x(j,:));
    f = [f; ff];
    FF = feval(Fun,x(j,:));
    if size(FF,1) > 1, FF = FF'; end
    F = [F; FF];
    fm(j,1) = softmerit(ff,FF,F1,F2,f0,Delta,sigma);
    if ff <= fglob & min(FF'-F1+sigma)>=0 & min(F2+sigma-FF')>=0
      ncall0 = ncall0 + j;
      xsoft = x(j,:);
      [fbestn,jbest] = min(fm(:,1));
      if fbestn < fbest
        fbest = fbestn;
        xbest = x(jbest,:);
      end
      ncall0,xbest,fbest,xsoft % show number of function values, best point
                               % and function value and the point xsoft
                               % fulfilling the conditions of the soft optimality
                               % theorem (hopefully close to xglob)
      return
    end
  end
  fm(:,2) = sqrt(eps)*ones(size(request,1),1); 
  ncall0 = ncall0 + size(request,1); % update function call counter
  [fbestn,jbest] = min(fm(:,1)); % best of the new merit function values
  % if a better function value has been found, update fbest
  if fbestn < fbest
    fbest = fbestn;
    xbest = x(jbest,:);
    ncall0,xbest,fbest % display current number of function values, best
                       % point and best function value if the latter
                       % have changed
  end
  if fbest < 0 & change == 0
    K = size(x1,1);
    ind = find(min(F-ones(K,1)*F1',[],2)>-eps&min(ones(K,1)*F2'-F,[],2)>-eps);
    if ~isempty(ind)
      change = 1;
      f0 = min(f(ind));
      Delta = median(abs(f-f0));
      disp('Change the merit function')
      for j=1:K
        fm(j,1) = softmerit(f(j),F(j,:),F1,F2,f0,Delta,sigma);
      end
      fm(:,2) = sqrt(eps)*ones(K,1);
      x = x1;
    end
  end    
end
ncall0,xbest,fbest  % show number of function values, best point and
                    % function value
