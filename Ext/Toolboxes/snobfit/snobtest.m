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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobtest.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file for running SNOBFIT on a set of test functions
% the test functions and their default box bounds are in the 
% subdirectory testfun
clear, clear mex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beginning of data to be adapted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = 'test';  % the intermediate data (after each call to SNOBFIT) 
                % are stored in test.mat
fcn = 'sh10';   % Shekel 10 function   
% fcn in {'bra','cam','gpr','hm3','hm6','ros','sh10','sh5','sh7',shu'}
fac = 0;        % factor for multiplicative perturbation of the data
ncall = 10000;   % limit on the number of function calls
[u,v,fglob] = defaults(fcn); % default box bounds, minimum of fcn
n = length(u);  % dimension of the problem
% the following are meaningful default values
npoint = n+6;   % number of random start points to be generated
nreq = n+6;     % no. of points to be generated in each call to SNOBFIT
x = rand(npoint,n);
x = x*diag(v-u) + ones(npoint,1)*u'; % generation of npoint random
                % starting points in [u,v]
dx = (v-u)'*1.e-5; % resolution vector
p = 0.5;        % probability of generating a point of class 4
prt = 0;        % print level
                % prt = 0 prints ncall, xbest and fbest if xbest has
                %         changed
                % prt = 1 in addition prints the points suggested by
                %         SNOBFIT, their model function values and 
                %         classes after each call to SNOBFIT        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of data to be adapted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:npoint
  f(j,:) = [feval(fcn,x(j,:))+fac*randn max(sqrt(eps),3*fac)];
% computation of the function values (if necessary, with additive
% noise)
end
ncall0 = npoint;   % function call counter
params = struct('bounds',{u,v},'nreq',nreq,'p',p); % input structure
% repeated calls to Snobfit
while ncall0 < ncall % repeat till ncall function values are reached
% (if the stopping criterion is not fulfilled first)
  if ncall0 == npoint  % initial call
    [request,xbest,fbest] = snobfit(file,x,f,params,dx);
    ncall0,xbest,fbest
  else                 % continuation call
    [request,xbest,fbest] = snobfit(file,x,f,params);
  end
  if prt>0, request, end
  clear x
  clear f
  for j=1:size(request,1)
    x(j,:) = request(j,1:n);
    f(j,:) = [feval(fcn,x(j,:))+fac*randn max(sqrt(eps),3*fac)];
  % computation of the (perturbed) function values at the suggested points
  end 
  ncall0 = ncall0 + size(f,1); % update function call counter
  [fbestn,jbest] = min(f(:,1)); % best function value
  if fbestn < fbest
    fbest = fbestn;
    xbest = x(jbest,:);
    ncall0,xbest,fbest % display current number of function values,
                       % best point and function value if fbest has
                       % changed
  end
  % check stopping criterion 
  % if fglob == 0, stop if fbest < 1.e-5
  % otherwise, stop if (fbest-fglob)/abs(fglob) < 1.e-2
  if fglob 
    if abs((fbest-fglob)/fglob) < 1.e-2,break,end
  else
    if abs(fbest) < 1.e-5,break,end
  end
end
ncall0,xbest,fbest  % show number of function values, best point and
                    % function value

