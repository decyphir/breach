function [xs, traj] = GetTrajSensi(Sys, traj, i_vars, i_params)
%GETTRAJSENSI gets the sensitivity of variables i_vars with respect to the
% parameters i_params.
% 
% Synopsis: [xs, traj] = GetTrajSensi(Sys, traj, i_vars, i_params)
% 
% Inputs:
%  - Sys      : the system.
%  - traj     : the trajectory based on which we compute the sensitivity.
%  - i_vars   : indexes or names of the variables. Unvalid variable names
%               or indexes will be skipped.
%  - i_params : indexes or names of the parameters. Parameter names not in
%               Sys.ParamList or parameter indexes not included in
%               [1, numel(traj.param)] will be skipped.
% 
% Outputs:
%  - xs   : array of dimension numel(i_vars) x numel(traj.time).
%  - traj : the initial field traj augmented of the computed sensitivity
% 
% Example (Lorentz84):
%   
%

% manage inputs
if ~isnumeric(i_vars)
    i_vars = FindParam(Sys,i_vars);
end
i_vars = i_vars(i_vars<=Sys.DimX);
i_vars = i_vars(i_vars>0);

if ~isnumeric(i_params)
    i_params = FindParam(Sys,i_params);
    i_params = i_params(i_params<=numel(Sys.p));
else
    i_params = i_params(i_params<=numel(traj.param));
end
i_params = i_params(i_params>0);

% check if it was already computed
if isfield(traj, 'sensis')
    sensis = traj.sensis;
else
    sensis = [];
    traj.XS = [];
end

iis = find(sensis==i_params); % check if the sensitivity is already computed
if ~isempty(iis)
    xs = traj.XS((iis-1)*Sys.DimX+i_vars,:);
    idx = traj.time>-1e90 & traj.time<1e90;
    xs = interp1(traj.time(idx), xs(idx), traj.time, 'linear', 0);  %NM: What's the point ??!!
    return
end

% recompute

% beurk.
time = traj.time(traj.time>-1e90);
time = time(time<1e90);

Ptmp = CreateParamSet(Sys, i_params);
Ptmp.pts = traj.param';
Ptmp = ComputeTrajSensi(Sys, Ptmp, time); % compute the sensitivity of variables to i_params

traj.sensis = [sensis i_params];
XS = Ptmp.traj{1}.XS;
XS = interp1(time, XS', traj.time, 'linear', 0);
if(Sys.DimX~=1) % if we have more than one variable, we have to transpose XS
    XS = XS';
end
traj.XS = [traj.XS; XS];  %NM: whould it be better to do: if any(traj.time>1e90) traj.XS = [traj.XS; [Ptmp.traj{1}.XS(:,1) Ptmp.traj{1}.XS Ptmp.traj{1}.XS(:,end)]]; end
xs = XS(i_vars,:);

end
