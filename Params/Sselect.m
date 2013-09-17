function Pn = Sselect(P, kn)
%SSELECT extract parameters from indices
%
% Synopsis:  Pn = Sselect(P [, kn])
%
% Inputs:
%  - P : the parameter set from which we extract a sub-parameter set
%
%  - kn: (Optional, defaut = find(P.selected~=0)) the indices of the set of
%        values to extract from P. If kn is not given and P does not have a
%        field selected, Pn is P. If kn is empty, Pn is P and a warning is
%        thrown.
%
% Output:
%  - Pn: parameter set containing the kn set of values of P. If the field
%        traj_ref wasn't defined in P, it is created in Pn.
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

% Manage input
if exist('kn','var')
    if isempty(kn)
        warning('Sselect:EmptyKn','The optional parameter ''kn'' is empty. Pn is defined as P.');
        Pn = P;
        return ;
    elseif(isequal(kn,1)&&(size(P.pts,2)==1))
        Pn = P;
        return;
    end
else
    try
        kn = find(P.selected~=0);
    catch %#ok<CTCH>
        Pn = P;
        return ;
    end
end

% Do the copy
Pn = [];


field_list_copy = {'dim', 'ParamList', 'DimX', 'DimP', 'props_names', 'props'};

for ii = 1:numel(field_list_copy)
    if isfield(P, field_list_copy{ii})
        Pn.(field_list_copy{ii}) = P.(field_list_copy{ii});
    end
end


field_list_select= {'pts', 'epsi', 'selected', 'props_values'};

for ii = 1:numel(field_list_select)
    if isfield(P, field_list_select{ii})
        Pn.(field_list_select{ii}) = P.(field_list_select{ii})(:,kn);
    end
end


if isfield(P,'traj')
    if ~isfield(P, 'traj_ref')
        P.traj_ref = 1:numel(P.traj);
    end
end

field_list_traj_ref_select = {'traj', 'etraj'};

for ii = 1:numel(field_list_traj_ref_select)
    if isfield(P, field_list_traj_ref_select{ii})
        Pn.(field_list_traj_ref_select{ii}) = P.(field_list_traj_ref_select{ii})(P.traj_ref(kn));
    end
end

if isfield(P, 'traj_to_compute')
    Pn.traj_to_compute = intersect(kn,P.traj_to_compute,'sorted');
else
    %TODO : to improve to avoid computation of two identical simulations
    Pn.traj_to_compute = 1:size(Pn.pts,2);
end

%TODO : set Pn.traj_ref

field_list_traj_ref_select2 = {'XS0', 'Xf', 'ExpaMax', 'XSf'};

for ii = 1:numel(field_list_traj_ref_select2)
    if isfield(P, field_list_traj_ref_select2{ii})
        Pn.(field_list_traj_ref_select2{ii}) = P.(field_list_traj_ref_select2{ii})(:,P.traj_ref(kn));
    end
end

end
