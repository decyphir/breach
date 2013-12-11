function P = CreateParamSet(Sys, Param, Ranges, Nb_pts)
%CREATEPARAMSET creates a parameter set structure for a system.
% 
% Synopsis: P = CreateParamSet(Sys[, Param[, Ranges[, Nb_pts]]])
% 
% Inputs:
%  -  Sys    : System under study
%  -  Param  : (Optional, default=all parameters) list of varying
%              parameters names or indexes. If it contains a parameter name
%              not existing in Sys.ParamList, this parameter is added to
%              the generated parameter set. The value of these new
%              parameters is either the mean value of the interval in the
%              Ranges is provided, zero otherwise.
%  -  Ranges : (Optional, default=[0.9*Sys.p, 1.1*Sys.p] for Sys.p~=0, 1
%              otherwise) array of size numel(Param) x 2 describing
%              parameter ranges.
%  -  Nb_pts : (Optional, default=1) number of parameter vectors in the
%              created parameter set;   same argument as Refine function
% 
% Outputs:
%  -  P : An initial parameter set with one or many parameter vectors
% 
% Example 1:
%     CreateSystem;
%     P = CreateParamSet(Sys); % x0, x1, x2, a, b, F and G are unknown
%                              % parameters. There values are (resp.) 0, 0,
%                              % 0, 0.25, 4, 0.5 and 0.5 and the epsi are
%                              % 1, 1, 1, 0.025, 0.1 and 0.05 .
% 
% Example 2:
%     CreateSystem;
%     P = CreateParamSet(Sys,'a'); % only a is an unknown parameter
% 
% Example 3:
%     CreateSystem;
%     P = CreateParamSet(Sys,{'a','x0'},[0.2,0.4 ; -0.6,0.6]); % a and x0
%                                      % are unknown parameters. The value
%                                      % of a % is 0.3, with an epsi of
%                                      % 0.1, and % the value of x0 is 0,
%                                      % with an epsi % of 0.6
% 
%See also: Refine SetParam SetEpsi SAddUncertainParam ComputeTraj
%

pts = Sys.p;

if exist('Param','var')
    if ischar(Param)
        Param = {Param};
    end
    nbParam = numel(Param);
    if isnumeric(Param)
        dim = reshape(Param,1,nbParam); % ensure that dim is a line vector
        ParamList = Sys.ParamList;
    else
        dim = FindParam(Sys,Param);
        ParamList = [Sys.ParamList{:}, Param(dim>size(pts,1))]; % extend ParamList to include formula parameters
    end
    %Here dim contains the indexes in pts of the varying parameters
    
    pts(dim(dim>size(pts,1))) = 0; %initialize formula parameters to 0
    
    epsi = zeros(nbParam,1);
    if exist('Ranges','var')
        if(size(Ranges,1)==2 && size(Ranges,2)==nbParam && nbParam~=2)
            warning('CreateParamSet:badRangesSize','Caution: size(Ranges) must be numel(Params) x 2. Fixed.');
            Ranges = Ranges'; % in case Ranges is in the wrong direction
        elseif(size(Ranges,2)~=2)
            error('CreateParamSet:wrongRangesSize','The size of the Ranges argument must be numel(Params) x 2.');
        end
        pts(dim) = (Ranges(:,2)+Ranges(:,1))/2;
        epsi(:,1) = Ranges(:,2)-pts(dim);
    else
        ptsun = pts(dim);
        epsi(ptsun~=0,1) = abs(ptsun(ptsun~=0)/10);
        epsi(ptsun==0,1) = 1;
    end
    
elseif(~isfield(Sys,'type')||(Sys.DimX==Sys.DimP))
    dim = 1:size(pts,1); %All the parameters are varying
    ParamList = Sys.ParamList;
    epsi = zeros(size(pts,1),1);
    epsi(pts~=0,1) = abs(pts(pts~=0)/10);
    epsi(pts==0,1) = 1;
    
else
    switch(Sys.type)
        case 'Breach'
            dim = 1:size(pts,1); %All the parameters are varying
        otherwise
            dim = Sys.DimX+1; %first parameter is varying
    end
    ptsun = pts(dim);
    epsi = zeros(numel(dim),1);
    epsi(ptsun~=0,1) = abs(ptsun(ptsun~=0)/10);
    epsi(ptsun==0,1) = 1;
    ParamList = Sys.ParamList;
end

P.DimX = Sys.DimX;
P.DimP = Sys.DimP;
P.dim = dim;
if(size(pts,2)>size(pts,1)) % should not happen !
    pts = pts'; % we want a column vector
end
if(size(epsi,2)>size(epsi,1)) % should not happen !
    epsi = epsi'; % we want a column vector
end

P.pts = pts;
P.epsi = epsi;
P.ParamList = ParamList;

P.traj_ref = 0;
P.traj_to_compute = 1;

if exist('Nb_pts', 'var') 
    P = Refine(P, Nb_pts);
end

end
