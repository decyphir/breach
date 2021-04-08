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
%%%%%%%%%%%%%%%%%%%%%%%%%% snobsoftdriver.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% driver for applying the soft optimality theorem and SNOBFIT to a soft
% constrained problem
clear; clear mex



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beginning of data to be adapted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% problem specification
file = 'hs';          % filename for storing intermediate data;
                      % after each call to SNOBFIT) are stored 
                      % in <file>.mat
fun = 'hsf18';        % name of objective function 
	              % of the form function f = fun(x) 
Fun = 'hsF18';        % name of constraint function
	              % of the form function F = Fun(x)
	              % here: Hock-Schittkowski problem no. 18
addpath('hsfun');     % add path to objective function
u = [2; 0];           % lower bounds for the variables
v = [50; 50];         % upper bounds for the variables
F1 = [0; 0];          % lower bounds for the constraints
F2 = [Inf; Inf];      % upper bounds for the constraints
sigma = 0.05*[25;25]; % softness degrees of constriants 
n = length(u);        % problem dimension (number of variables)
m = length(F1);       % number of constraints 

% meaningful default values for SNOBFIT parameters
nreq = n+6;           % number of points to be generated in each call 
                      % to SNOBFIT
dx = (v-u)'*1.e-5;    % resolution vector
p = 0.5;              % probability of generating a point of class 4
prt = 0;              % print level
                      % prt = 0 prints ncall, xbest and fbest 
                      % if xbest has changed
                      % prt = 1 in addition prints the points suggested
                      % by SNOBFIT, their model function values and 
                      % classes after each call to SNOBFIT  

% parameters used in this driver      
% (The default for minfcall should be chosen such that the the local 
% quadratic fit around the best point uses at least once the maximal 
% number of points, which requires minfcall >= n*(n+3)+1+nreq.
% In cases where function evaluations are not very expensive, 
% minfcall should be chosen much higher.)
ncall = 10000;        % limit on the number of function calls
minfcall = 500;       % minimum number of function values before
                      % considering stopping
nstop = 5;            % number of times no improvement is tolerated
% the optimization is stopped either if approximately ncall fucntion 
% values have been exceeded, or if at least minfval function values 
% were obtained and the best function value wasn't improved in the last
% nstop calls to SNOBFIT  
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of data to be adapted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



params = struct('bounds',{u,v},'nreq',nreq,'p',p); % input structure
[request,xbest,fbest] = snobfit(file,[],[],params,dx);
% initial call with empty list of input points
x1 = [];
f = [];
F = [];
change = 0;
% computation of the objective function values and constraints for the 
% points suggested by SNOBFIT
f0 = Inf; % initialization of f0
npoint=size(request,1);
for j=1:npoint
  x(j,:) = request(j,1:n);
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
% computation of the merit function for the points
for j=1:npoint
  fm(j,1) = softmerit(f(j),F(j,:),F1,F2,f0,Delta,sigma);
end
fm(:,2) = sqrt(eps)*ones(npoint,1); % set uncertainty of the merit 
                                    % function
ncall0 = size(request,1); % function call counter
[fbest,jbest] = min(fm(:,1)) % best of the new merit function values]
xbest = x(jbest,:);
ncall0,xbest,fbest % display current number of function values, best
                   % point and best merit function value
nstop0 = 0;
while ncall0 < ncall  % repeat until the limit on function calls is 
                      % reached
  [request,xbest,fbest] = snobfit(file,x,fm,params);
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
    nstop0 = 0;
  elseif ncall >= minfcall,
    nstop0 = nstop0 + 1;
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
  % check stopping criterion
  if nstop0 >= nstop & ncall0 >= minfcall, break, end 
end
ncall0,xbest,fbest  % show number of function values, best point and
                    % function value
softfeas = [];
fmerit = []'
for j=1:size(F,1)
  delta(j,:) = max(0,max(F1'-F(j,:),F(j,:)-F2')./sigma');
  if max(delta(j,:))<=1 
    softfeas = [softfeas; j];
    fmerit0 = softmerit(f(j),F(j,:),F1,F2,f0,Delta,sigma);
    fmerit = [fmerit; fmerit0];
  end
end
if ~isempty(fmerit)
  [fsoft,isoft] = min(fmerit);
  xsoft = x1(softfeas(isoft),:); 
  ncall0,xbest,fbest,xsoft
else
  delta0 = sum(delta.^2,2);
  [mindelta,isoft] = min(delta0);
  xsoft = x1(isoft,:);
  ncall0,xbest,fbest,xsoft
  disp('WARNING: No soft feasible point has been found')
  disp('xsoft is the point with smallest constraint violation sum(delta_i^2')
end	



