function Pn = Sselect(P, kn)
%SSELECT extracts parameters from indices
% 
% Synopsis:  Pn = Sselect(P[, kn])
% 
% Inputs:
%  - P  : the parameter set from which we extract a sub-parameter set
%  - kn : (Optional, defaut=find(P.selected~=0)) the indices of the
%         parameter vectors to extract from P. If kn is not given and P
%         does not have a field selected, Pn is P. If kn is empty, Pn is P
%         and a warning is thrown.
% 
% Output:
%  - Pn : parameter set containing the selected required sub-parameter set.
%         If the fields traj_ref or traj_to_compute wasn't defined in P,
%         they are created in Pn.
% 
% Example (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys, {'a','F'}, [0.15,0.35;5,25], 4); % creates 16 param vectors
%   
%   Peven = Sselect(P, 2:2:16);
%   Podd = Sselect(P, 1:2:15);  % split P into two param sets
%   
%   P = ComputeTraj(Sys, P, 0:0.01:10);
%   phi = STL_Formula('phi', 'ev (x1[t]<-3)');
%   [~,val] = SEvalProp(Sys,P,phi,0);
%   Pvalid = Sselect(P,find(val>0));  % the four select param vectors verifying phi
%   figure ; SplotVar(Pvalid)
%   
%   P2 = CreateParamSet(Sys,'a',[0,2]);
%   P2 = ComputeTraj(Sys,P2,1:0.1:10);
%   P2 = SetParam(P2,'propParam',2);
%   P2 = SAddUncertainParam(P2,'propParam');
%   P2 = Refine(P2,[3,3]);
%   P2.traj_ref   % have three parameter vector leading to the trajectory
%   P2 = Sselect(P2,[2,5]);  % select two of them
%   size(P2.traj,2) % only one traj
%   P2.traj_ref     % with both parameter vector leading to it
%   P2.traj_to_compute  % and no parameter vector remaining to compute
% 
%See also SConcat PInclude CreateParamSet SAddUncertainParam
%

% val = STL_Eval(Sys, phi, P, P.traj, 0);     % this code also works for
% Pvalid = Sselect(P, P.traj_ref(find(val>0))); % selecting valid parameter sets

% Manage input
if exist('kn','var')
    kn = kn(kn>0); % keep only valid kn
    kn = kn(kn<=size(P.pts,2));
    if isempty(kn)
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
        return;
    end
end

% Do the copy
Pn = struct();


field_list_copy = {'dim', 'ParamList', 'DimX', 'DimP', 'props_names', 'props', 'time_mult'};

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


field_list_traj_ref_select = {'traj', 'etraj'};
field_list_traj_ref_select2 = {'XS0', 'Xf', 'ExpaMax', 'XSf'};

Pn.traj_ref = zeros(1,numel(kn));

if ~isfield(P, 'traj_ref')
    for ii = 1:numel(field_list_traj_ref_select)
        if isfield(P, field_list_traj_ref_select{ii})
            Pn.(field_list_traj_ref_select{ii}) = P.(field_list_traj_ref_select{ii})(kn);
        end
        if isfield(P, field_list_traj_ref_select2{ii})
            Pn.(field_list_traj_ref_select2{ii}) = P.(field_list_traj_ref_select2{ii})(:,kn);
        end
    end
    if isfield(Pn,'traj')
        idx_max = min(numel(Pn.traj),numel(kn));
        Pn.traj_ref(1:idx_max) = 1:idx_max;
    end
else
    kn_ref = P.traj_ref(kn);
    kn_ref = unique(kn_ref(kn_ref~=0),'stable'); % don't copy 1/ non-computed traj OR 2/ many times same stuff
    for ii = 1:numel(field_list_traj_ref_select)
        if isfield(P, field_list_traj_ref_select{ii})
            Pn.(field_list_traj_ref_select{ii}) = P.(field_list_traj_ref_select{ii})(kn_ref);
        end
        if isfield(P, field_list_traj_ref_select2{ii})
            Pn.(field_list_traj_ref_select2{ii}) = P.(field_list_traj_ref_select2{ii})(:,kn_ref);
        end
    end
    for ii=1:numel(kn_ref) % to set the correct values to traj_ref, we do a mapping from the unique...
        Pn.traj_ref(P.traj_ref(kn)==kn_ref(ii)) = ii; % values in kn_ref to 1:numel(unique(kn_ref))
    end
end

[~,Pn.traj_to_compute] = unique(Pn.pts(1:Pn.DimP,:)','rows','first'); % set traj_to_compute
if isfield(Pn,'traj') % optimizing test
    Pn.traj_to_compute = setdiff(Pn.traj_to_compute,find(Pn.traj_ref~=0)); % don't keep those already computed
end
Pn.traj_to_compute = sort(reshape(Pn.traj_to_compute,1,[])); % set it in a line shape

end
