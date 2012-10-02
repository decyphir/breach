function S = CreateParamSet(Sys,Param,Ranges)
%   CREATEPARAMSET Create a parameter set structure for a system. If the
%   argument Ranges is not provided, the range is equal to 1/10th of the
%   parameter value, except if the parameter value is 0, in which case, the
%   range is set to 1.
%
%   Synopsis: P = CreateParamSet(Sys [,Param,Ranges] )
%
%   Inputs:
%
%    -  Sys      System under study
%    -  Param    List of varying parameters
%    -  Ranges   Parameters ranges
%
%   Outputs:
%
%    -  P        An initial parameter set with one point in the center of
%                ranges of the varying paramaters and the default value
%                in Sys.p
%

pts = Sys.p;

if (exist('Param','var'))
    if ischar(Param)
        Param = {Param};
    end
    nbParam = numel(Param);
    %NM dim = zeros(nbParam,1);
    if isnumeric(Param)
        dim = Param(1:nbParam); % ensure that dim is a vector
        ParamList = Sys.ParamList;
    else
        ind = FindParam(Sys,Param);
        new_params = Param(ind>Sys.DimP);
        ParamList = [Sys.ParamList{:} new_params];
        dim = ind;
    end
    
    pts(dim)=0; %initialize all parameters (new and not new) values to 0
    
    dim_sys = dim(dim<=Sys.DimP);
    pts(dim_sys) = Sys.p(dim_sys); % copy the not new parameters
    
    ptsun = pts(dim);
    epsi = zeros(nbParam,1);
    if (exist('Ranges','var'))
        pts(dim) = (Ranges(:,2)+Ranges(:,1))/2;
        ptsun = pts(dim);
        epsi(:,1) = Ranges(:,2)-ptsun;
    else
        epsi(ptsun~=0) = abs(ptsun(ptsun~=0)/10);
        epsi(ptsun==0) = 1;
    end
    
else
    dim = 1:Sys.DimP;
    ptsun = pts(dim);
    epsi = zeros(numel(dim),1);
    epsi(ptsun~=0) = abs(ptsun(ptsun~=0)/10);
    epsi(ptsun==0) = 1;
    ParamList = Sys.ParamList;
end

S.DimX = Sys.DimX;
S.DimP = Sys.DimP;
S.dim = dim;
if (size(pts,2)>size(pts,1))
    pts = pts';
end
if (size(epsi,2)>size(epsi,1))
    epsi = epsi';
end

S.pts = pts;
S.epsi = epsi;
S.ParamList=ParamList;
end
