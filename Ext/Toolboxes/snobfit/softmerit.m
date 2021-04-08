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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% softmerit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fm = softmerit(f,F,F1,F2,f0,Delta,sigma)
% merit function of the soft optimality theorem
%
% Input:
% f        objective function value
% F        vector containing the values of the constraint functions
%          (m-vector)
% F1       m-vector of lower bounds of the constraints
% F2       m-vector of upper bounds of the constraints
% f0       scalar parameter in the merit function
% Delta    scalar, positive parameter in the merit function
% sigma    positive m-vector, where sigma(i) is the permitted violation 
%          of constraint i
%
function fm = softmerit(f,F,F1,F2,f0,Delta,sigma)
if ~isfinite(f) | any(~isfinite(F)) 
% if the objective function or one of the constraint functions is 
% infinite or NaN, set the merit function value to 3
  fm = 3;
  return
end
m = length(F);
delta = 0;
for i=1:m
  if F(i) < F1(i)
    delta = delta + (F1(i)-F(i))^2/sigma(i)^2;
  elseif F(i) > F2(i)
    delta = delta + (F(i)-F2(i))^2/sigma(i)^2;
  end
end
fm = (f-f0)/(Delta+abs(f-f0)) + 2*delta/(1+delta);
