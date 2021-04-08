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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snob5.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x1 = snob5(x,u,v,dx,nreq)
% generates nreq points of class 5 in [u,v]
%
% Input:
% x       the rows are the points already chosen by Snobfit (can be 
%         empty)
% u,v     bounds of the box in which the points are to be generated
% dx      resolution vector, i.e. the ith coordinate of a point to be
%         generated is an integer-valued multiple of dx(i) 
% nreq    number of points to be generated
%
% Output:
% x1      the rows are the requested points
%         x1 is of dimension nreq x n (n = dimension of the problem)



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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobinput.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,np,t] = snobinput(x,f)
% checks whether there are any duplicates among the points given by the
% rows of x, throws away the duplicates and computes their average
% function values and an estimated uncertainty
%
% Input:
% x	the rows of x are a set of points
% f	f(j,1) is the function value of x(j,:) and f(j,2) is its
%	uncertainty
%
% Output:
% x	updated version of x (possibly some points have been deleted)
% f	updated version of f (f(j,1) is the average function value and
%	f(j,2) the estimated uncertainty pertaining to x(j,:))  
% np	np(j) is the number of times the row x(j,:) appeared in the
%	input version of x
% t	t(j) is np(j) times the variance of the function values measured
%	for point x(j,:)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snoblocf.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [y,f,g,sigma] = snoblocf(j,x,f,near,dx,u,v)
% computes a local fit around the point x0 = x(j,:) and minimizes it
% on a trust region
%
% Input:
% j         index of the point around which the fit is to be computed 
% x         the rows contain the points where the function has been
%           evaluated
% f         the corresponding function values and their uncertainties,
%           i.e., f(j,1) = f(x(j,:)), f(j,2) = df(x(j,:))
% near      near(j,:) is a vector containing the indices of the nearest
%           neighbors of the point x(j,:)
% dx        resolution vector, i.e. the ith coordinate of a point to be
%           generated is an integer-valued multiple of dx(i) 
% u,v       bounds of the box where the points should be generated
%
% Output:
% y         estimated minimizer in the trust region
%           
% f1        its estimated function value
% g         estimated gradient for the fit
% sigma     sigma = norm(A*g-b)/sqrt(K-n), where A and b are the 
%           coefficients resp. right hand side of the fit, n is the
%           dimension and K the number of nearest neighbors considered
%           (estimated standard deviation of the model errors)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snoblp.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [local,nlocal] = snoblp(f,near,ind)
% computes a pointer to all `local' points (i.e. points whose neighbors
% have `significantly larger' function values)
%
% Input:
% f	f(j) is the function value of point j
% near	near(j,:) contains the indices of the size(near,2) neighbors of
%	point j
% ind	pointer to the boxes to be considered (optional, default
%       1:length(f))
%
% Output:
% local		vector containing the indices of all local points
% nlocal	vector containing the indices of all nonlocal points



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobnan.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = snobnan(fnan,f,near,inew)
% replaces the function values NaN of a set of points by a value 
% determined by their nearest neighbors with finite function values, 
% with a safeguard for the case that all neighbors have function value
% NaN
%
% Input:
% fnan	  	vector containing the pointers to the points where the
%               function value could not be obtained
% f	  	f(:,1) set of available function values
%		f(:,2) their uncertainty/variation
% near(j,:)	vector pointing to the nearest neighbors of point j
% inew          vector pointing to the new boxes and boxes whose nearest
%               neighbors have changed
%
% Output:
% f		updated version of f



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobnewb.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xl,xu,vol,u,v] = snobnewb(xnew,xl,xu,small,u,v,u1,v1)
% if xnew contains points outside the box [u,v], the box is made larger
% such that it just contains all new points; moreover, if necessary, 
% [u,v] is made larger to contain [u1,v1]
% the box bounds xl and xu of the boxes on the boundary and the volume 
% measures small of these boxes are updated if necessary
%
% Input:
% xnew	the rows of xnew contain the new points
% xl	xl(j,:) is the lower bound of box j
% xu	xu(j,:) is the upper bound of box j
% small	small(j) is the volume measure of box j
% u,v	old box bounds
% u1,v1 box in which the points are to be generated
%
% Output:
% xl	updated version of xl
% xu 	updated version of xu
% small	updated version of small
% u,v	updated box bounds (the old box is contained in the new one)



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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobpoint.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = snobpoint(x,xl,xu,f0,g,sigma,u,v,dx)
% for a box [xl,xu] containing a point x, a point y in the intersection 
% of [xl,xu] and [u,v] is constructed such that it is both not close to 
% x and to the boundary of [xl,xu] and its function value is estimated
% from a local quadratic model around x
%
% Input:
% x	     point contained in [xl,xu]
% xl,xu	     box bounds
% u,v        the point is to be generated in [u,v]
% f0         f0(1) is the function value at x, f0(2) is its uncertainty
% g,G,sigma  the local quadratic model around x is given by
%            q(y)=f0(1)+g*(y-x)'+sigma*((y-x)*diag(D)*(y-x)'+f0(2))
%            for a row vector y, where D = f0(2)./dx.^2
% dx         resolution vector
%
% Output:
% y	     point in the intersection of [xl,xu] and [u,v]
% f          corresponding estimated function value



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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobqmin.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x = snobqmin(a,b,xl,xu)
% minimization of the quadratic polynomial p(x) = a*x^2+b*x over [xl,xu]
%
% Input:
% a, b    coefficients of the polynomial
% xl,xu   bounds (xl < xu)
% 
% Output:
% x       minimizer in [xl,xu]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobround.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x = snobround(x,u,v,dx)
% a point x is projected into the interior of [u,v] and x(i) is
% rounded to the nearest integer multiple of dx(i)
%
% Input:
% x	vector of length n
% u,v	vectors of length n such that u < v
% dx    vector of length n
%
% Output:
% x	projected and rounded version of x



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% snobsoftdriver.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% driver for applying the soft optimality theorem and SNOBFIT to a soft
% constrained problem



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% snobsofttest.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file (m-script) for applying the soft optimality theorem and SNOBFIT
% to some Hock-Schittkowski problems



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobsplit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xl,xu,x,f,nsplit,small] = snobsplit(x,f,xl0,xu0,nspl,u,v)
% splits a box [xl0,xu0] contained in a bigger box [u,v] such that
% each of the resulting boxes contains just one point of a given set of
% points
%
% Input:
% x	the rows are a set of points
% f	f(j,:) contains the function value, its variation and possibly
%	other parameters belonging to point x(j,:)
% xl0	vector of lower bounds of the box
% xu0	vector of upper bounds of the box
% nspl	nspl(i) is the number of splits the box [xl0,xu0] has already
%	undergone along the ith coordinate
% u,v	bounds of the original problem	
%	[xl0,xu0] is contained in [u,v]
%
% Output:
% xl	xl(j,:) is the lower bound of box j
% xu	xi(j,:) is the upper bound of box j
% x      	x(j,:) is the point contained in box j
% f	f(j,:) contains the function value at x(j,:), its uncertainty
%	etc. 
% nsplit nsplit(j,i) is the number of times box j has been split in the
%	ith coordinate
% small	small(j) = integer-valued logarithmic volume measure of box j



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobtest.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file for running SNOBFIT on a set of test functions
% the test functions and their default box bounds are in the 
% subdirectory testfun



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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% snobwarn.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% issues a warning if SNOBFIT has not been able to generate the desired
% number of points



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




