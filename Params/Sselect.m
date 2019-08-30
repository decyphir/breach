function P = Sselect(P, kn)
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
        return ;
    elseif(isequal(kn,1)&&(size(P.pts,2)==1))
        return;
    end
else
    try
        kn = find(P.selected~=0);
    catch %#ok<CTCH>
        return;
    end
end

field_list_select= {'pts', 'epsi', 'selected', 'props_values'};

for ii = 1:numel(field_list_select)
    if isfield(P, field_list_select{ii})
        P.(field_list_select{ii}) = P.(field_list_select{ii})(:,kn);
    end
end

%% Optimize for one
if numel(kn)==1
    if isfield(P, 'traj')
        if numel(P.traj)>=P.traj_ref(kn) % traj exist for this guy
            P.traj= P.traj(P.traj_ref(kn));
            P.traj_ref = 1;
            P.traj_to_compute = [];
        else
            P = rmfield(P, 'traj');
            P.traj_ref = 1;
            P.traj_to_compute = 1;
        end
    else
       P.traj_ref = 1;
       P.traj_to_compute = 1; 
    end
else
    P = Preset_traj_ref(P);
end

end
