function P = SCreateNull(P)
%SCREATENULL creates an empty parameter set
%
% Synopsis: P = SCreateNull([P])
%
% Input:
%  - P : a parameter set. If not provided, SCreateNull answers [].
%
% Output:
%  - P : the parameter set with fields pts, traj, epsi, ToCheck, XS0, Xf,
%        ExpaMax ans XSf set to []
%
%See also CreateParamSet
%

if(nargin==0)
    P=[];
end

P.pts = [];
P.traj = [];
P.epsi = [];
P.ToCheck = [];
P.XS0 = [];
P.Xf = [];
P.ExpaMax = [];
P.XSf = [];

end
