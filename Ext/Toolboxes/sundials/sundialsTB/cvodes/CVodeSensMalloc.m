function CVodeSensMalloc(Ns, meth, yS0, varargin)
%CVODESENSMALLOC allocates and initializes memory for FSA with CVODES.
% 
% Synopsis: CVodeSensMalloc(Ns, meth, yS0[, options])
% 
% Inputs:
%  - Ns      : is the number of parameters with respect to which
%              sensitivities are desired
%  - meth    : FSA solution method [ 'Simultaneous' | {'Staggered'} ]
%              Specifies the FSA method for treating the nonlinear system
%              solution for sensitivity variables. In the simultaneous
%              case, the nonlinear systems for states and all sensitivities
%              are solved simultaneously. In the Staggered case, the
%              nonlinear system for states is solved first and then the
%              nonlinear systems for all sensitivities are solved at the
%              same time.
%  - yS0     : Initial conditions for sensitivity variables.
%              yS0 must be a matrix with N rows and Ns columns, where N is
%              the problem dimension and Ns the number of sensitivity
%              systems.
%  - options : is an (optional) set of FSA options, created with the
%              CVodeSetFSAOptions function.
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2009-06-05 16:26:17 $

mode = 2;

if(nargin < 3)
    error('CVodeSensMalloc:tooFewParameters','CVodeSensMalloc got less than three parameters.');
end

options = [];
if(nargin > 3)
    options = varargin{1};
end

cvm(mode,Ns,meth,yS0,options);

end
