function Pn = Sselect(P, kn)
%SSELECT extract parameters from indices
%
% Synopsis:  Pn = Sselect(P [, kn])
%
% Inputs:
%  - P : the parameter set from which we extract a sub-parameter set
%  - kn: (Optional) the indices of the set of values to extract from P. If
%         kn is not given, looks for a 'selected' field in P and returns
%         Sselect(P,find(P.selected~=0)) if this field exists or P itself
%         otherwize.
%
% Output:
%  - Pn: parameter set containing the kn set of values of P
%
% Example (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys, {'a','F'}, [0.15,0.35;5,25], 4); % creates 16 param sets
%
%   Peven = Sselect(P, 2:2:16);
%   Podd = Sselect(P, 1:2:15);  % split P into two param sets
%
%   P = ComputeTraj(Sys, P, 0:0.01:10);
%   phi = QMITL_Formula('phi', 'ev (x1[t]<-3)');
%   [~,val] = SEvalProp(Sys,P,phi,0);
%   Pvalid = Sselect(P,find(val>0));  % select param sets verifying phi
%

% val = QMITL_Eval(Sys, phi, P, P.traj, 0);     % this code also works for
% Pvalid = Sselect(P, P.traj_ref(find(val>0))); % selecting valid parameter sets

if(nargin==1)
    try
        kn = find(P.selected~=0);
    catch %#ok<CTCH>
        Pn = P;
        return;
    end
end

Pn = [];

field_list_copy = {'dim', 'ParamList', 'DimX', 'DimP', 'props_names', 'props'};

for i = 1:numel(field_list_copy)
    if isfield(P, field_list_copy{i})
        Pn.(field_list_copy{i}) = P.(field_list_copy{i});
    end
end

field_list_select= {'pts', 'epsi', 'selected', 'props_values'};

for i = 1:numel(field_list_select)
    if isfield(P, field_list_select{i})
        Pn.(field_list_select{i}) = P.(field_list_select{i})(:,kn);
    end
end

field_list_traj_ref_select = {'traj', 'etraj'};

if isfield(P,'traj')
    if ~isfield(P, 'traj_ref')
        P.traj_ref = 1:numel(P.traj);
    end
end

for i = 1:numel(field_list_traj_ref_select)
    if isfield(P, field_list_traj_ref_select{i})
        Pn.(field_list_traj_ref_select{i}) = P.(field_list_traj_ref_select{i})(P.traj_ref(kn));
    end
end

field_list_traj_ref_select2 = {'XS0', 'Xf', 'ExpaMax', 'XSf'};

for i = 1:numel(field_list_traj_ref_select2)
    if isfield(P, field_list_traj_ref_select2{i})
        Pn.(field_list_traj_ref_select2{i}) = P.(field_list_traj_ref_select2{i})(:,P.traj_ref(kn));
    end
end
