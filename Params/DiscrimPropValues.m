function Sets = DiscrimPropValues(P)
%DISCRIMPROPVALUES separates trajectory set P according to properties
% checked for P.
% 
% Synopsis: Sets = DiscrimPropValues(P)
% 
% Input:
%  - P a parameter set. It may contain many parameter vectors (actually,
%      the function is useless if there is only one parameter vector).
% 
% Output:
%  - Sets an array of parameter set. All parameter vector of a parameter
%         set of Sets has the same boolean evalution at the first time
%         instant at which the formula is evaluated for all computed
%         properties. The parameter sets are augmented by a field status.
%         If no properties has been checked, Sets if equal to P.
% 
%See also SEvalProp
%

if ~isfield(P,'props_values')
    Sets = P;
    return;
end

pr = P.props_values;
nb_prop = size(pr,1);

traj_status = zeros(1, size(pr,2));
for ii=1:nb_prop
    pr_ii = cat(1, pr(ii,:).val);
    pr_ii = (pr_ii(:,1)>0)'; % keep the boolean satisfaction at first evaluated time instant and set it in a line array
    traj_status = traj_status+pr_ii*2.^ii; % this is a recap of evaluation of all properties
end

status = traj_status;

Sets = [];
P.status = 0;
while ~isempty(status)
    status_current = status(1);
    statusi = find(traj_status==status_current); % get all traj with the same status
    status = status(status~=status_current); % remove all status_current to status
    Pn = Sselect(P, statusi);
    Pn.status = status_current;
    Sets = [Sets, Pn]; %#ok<AGROW>
end

end
