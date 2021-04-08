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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobdriver.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% driver for running SNOBFIT 
clear, clear mex



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beginning of data to be adapted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% problem specification
file = 'test';        % filename for storing intermediate data;
                      % after each call to SNOBFIT) are stored 
                      % in <file>.mat
fcn = 'sh10';         % objective function (here: Shekel 10)
                      % function f = fcn(x)    
addpath('testfun');   % add path to objective function
u = zeros(4,1);       % lower box bounds
v = 10*ones(4,1);     % upper box bounds
n = length(u);        % problem dimension (number of variables)

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
for j=1:size(request,1)
  x(j,:) = request(j,1:n);
  f(j,:) = [feval(fcn,x(j,:)) sqrt(eps)];
% computation of the function values at the suggested points and setting
% uncertainties 
end
ncall0 = size(f,1); % function call counter
[fbestn,jbest] = min(f(:,1)); % best function value
xbest = x(jbest,:);
ncall0,xbest,fbest % display current number of function values, best 
                   % point and function value
nstop0 = 0;
% repeated calls to Snobfit
while ncall0 < ncall % repeat till ncall function values are reached
% (if the stopping criterion is not fulfilled first)
  [request,xbest,fbest] = snobfit(file,x,f,params);
  if prt>0, request, end
  clear x
  clear f
  for j=1:size(request,1)
    x(j,:) = request(j,1:n);
    f(j,:) = [feval(fcn,x(j,:)) sqrt(eps)];
  % computation of the function values at the suggested points
  end 
  ncall0 = ncall0 + size(f,1); % update function call counter
  [fbestn,jbest] = min(f(:,1)); % best function value
  if fbestn < fbest
    fbest = fbestn;
    xbest = x(jbest,:);
    ncall0,xbest,fbest % display current number of function values,
                       % best point and function value if fbest has
                       % changed
    nstop0 = 0;
  elseif ncall >= minfcall,
    nstop0 = nstop0 + 1;
  end
  % check stopping criterion 
  if nstop0 >= nstop & ncall0 >= minfcall, break, end 
end
ncall0,xbest,fbest  % show number of function values, best point and
                    % function value

