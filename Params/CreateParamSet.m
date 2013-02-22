function P = CreateParamSet(Sys, Param, Ranges, Nb_pts)
%   CREATEPARAMSET Create a parameter set structure for a system. If the
%   argument Param contains parameters name not existing in the system,
%   they are added to the parameter set. The value of these parameters is
%   either the median point of the ranges, if the range if provided, or
%   zero.
%   If the argument Param is not provided, all the parameters are
%   considered as varying parameters.
%   If the argument Ranges is not provided, the epsi is equal to 1/10th of
%   the parameter value, except if the parameter value is 0, in which case,
%   the epsi is set to 1.
%
%   Synopsis: P = CreateParamSet(Sys [,Param,Ranges,Nb_pts] )
%
%   Inputs:
%
%    -  Sys      System under study
%    -  Param    List of varying parameters
%    -  Ranges   Parameters range
%    -  Nb_pts   Number of points in the parameter set (default: 1) in each
%                dimension;   same argument as Refine function
%               
%
%   Outputs:
%
%    -  P        An initial parameter set with one point in the center of
%                ranges of the varying parameters and the default value
%                in Sys.p
%   Example 1:
%     CreateSystem;
%     P = CreateParamSet(Sys); % x0, x1, x2, a, b, F and G are unknown
%                              % parameters. There values are (resp.) 0, 0,
%                              % 0, 0.25, 4, 0.5 and 0.5 and the epsi are
%                              % 1, 1, 1, 0.025, 0.1 and 0.05 .
%
%   Example 2:
%     CreateSystem;
%     P = CreateParamSet(Sys,'a'); % only a is an unknown parameter
%
%   Example 3:
%     CreateSystem;
%     P = CreateParamSet(Sys,{'a','x0'},[0.2,0.4;-0.6,0.6]); % a and x0 are
%                                      % unknown parameters. The value of a
%                                      % is 0.3, with an epsi of 0.1, and
%                                      % the value of x0 is 0, with an epsi
%                                      % of 0.6
%
%See also: Refine SetParam
%

pts = Sys.p;

if (exist('Param','var'))
    if ischar(Param)
        Param = {Param};
    end
    nbParam = numel(Param);
    if isnumeric(Param)
        dim = Param(1:nbParam); % ensure that dim is a vector (NM:even for Param(1:1)??)
        ParamList = Sys.ParamList;
    else
        ind = FindParam(Sys,Param);
        new_params = Param(ind>size(pts,1)); %in case param contains parameters not
                                             %in Sys.ParamList (aka:we create parameters)
        ParamList = [Sys.ParamList{:} new_params];
        dim = ind;
    end
    %Here dim contains the indexes in pts of the varying parameters
    
    pts(dim)=0; %initialize all parameters (new and not new) values to 0
    
    dim_sys = dim(dim<=size(Sys.p,1));
    pts(dim_sys) = Sys.p(dim_sys); % copy the not new parameters
    
    epsi = zeros(nbParam,1);
    if (exist('Ranges','var'))
        pts(dim) = (Ranges(:,2)+Ranges(:,1))/2;
        epsi(:,1) = Ranges(:,2)-pts(dim);
    else
        ptsun = pts(dim);
        epsi(ptsun~=0) = abs(ptsun(ptsun~=0)/10);
        epsi(ptsun==0) = 1;
    end
    
else
    dim = 1:size(pts,1); %All the parameters are varying
    ptsun = pts(dim);
    epsi = zeros(numel(dim),1);
    epsi(ptsun~=0) = abs(ptsun(ptsun~=0)/10);
    epsi(ptsun==0) = 1;
    ParamList = Sys.ParamList;
end

P.DimX = Sys.DimX;
P.DimP = Sys.DimP;
P.dim = dim;
if (size(pts,2)>size(pts,1))
    pts = pts';
end
if (size(epsi,2)>size(epsi,1))
    epsi = epsi';
end

P.pts = pts;
P.epsi = epsi;
P.ParamList=ParamList;

if exist('Nb_pts', 'var') 
    P = Refine(P, Nb_pts);
end