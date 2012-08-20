function S = CreateParamSet(Sys,Param,Ranges)
%   CREATEPARAMSET Create a parameter set structure for a system
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
  dim = zeros(nbParam,1);
  if isnumeric(Param)
    dim = Param(1:nbParam);
    ParamList = Sys.ParamList;
  else
    ind = FindParam(Sys,Param);
    ParamList = unique({Sys.ParamList{:}, Param{:}});         
    dim = ind;
  end
  
  pts(dim)=0;
  
  dim_sys = dim(dim<= Sys.DimP);
  pts(dim_sys)= Sys.p(dim_sys);
  
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
