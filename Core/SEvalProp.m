function [P, val] = SEvalProp(Sys, P, phis, taus, ipts, break_level, method, VERBOSE)
%SEVALPROP Eval property for previously computed trajectories
% 
% Synopsis: [P, val] = SEvalProp(Sys, P, phis[ , taus[, ipts[, break_level[, method[, VERBOSE]]]]])
% 
% Inputs:
%  - Sys         : System structure
%  - P           : Parameter set. It may contain many parameter vectors.
%                  All trajectories must be computed or an error is thrown.
%  - phis        : an array of STL property(ies)
%  - taus        : (Optional) Time point(s) when to estimate properties. If
%                  not provided, the formulas are evaluated at the first
%                  time point of the trajectories. It is a cell array of
%                  size 1 x numel(phis), the ith cell describing the time
%                  point(s) where the ith formula should be evaluated. Taus
%                  may also be a scalar, in which case, all the formula are
%                  evaluated at this time point. The third possiblity is
%                  that taus is an array. In this case: 1/ numel(phis) is
%                  equal to 1, so the formula is evaluated at all time
%                  points provided in taus ; 2/ numel(taus)==numel(phis),
%                  in which case, the ith formula is evaluated at the ith
%                  time point provided in taus ; 3/ otherwise, all formula
%                  are evaluated at all time points provided in taus. As
%                  the two last possibilities may be confusing, it is
%                  recommanded to use a cell array.
%  - ipts        : (Optional, default or empty=all parameter sets) Indices
%                  of parameter vectors for which the formulas are evaluated.
%  - break_level : (Optional, default=0) defines the depth of breaking of
%                  the formulas. If lower or equal to 1, it is ignored. If
%                  greater or equal to two, SEvalProp provides the
%                  evaluation of formulas in phis and all sub-formulas
%                  until the depth provided.
%  - method      : (Optional, default or empty='thom') string indicating
%                  the method which must be used to evaluate the formulas.
%                  It must be 'classic' or 'thom'.
%  - VERBOSE     : (Optional, default=false) boolean indicating if progress
%                  bar is shown.
% 
% Outputs:
%  - P   : param set with prop_names, prop and prop_values fields
%  - val : an array containing the quantitative satisfaction of properties
%          for each trajectory at the first time point of tau. The
%          dimension of val is numel(props) x numel(ipts)
% 
% Example (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys, 'F', [10, 20]);
%   P = SetParam(P, 'the_threshold', 2);
%   P = ComputeTraj(Sys, P, 0:.01:10);
%   phi = STL_Formula('phi','ev_[0,1] (x0[t]>the_threshold)');
%   [~,val] = SEvalProp(Sys, P, phi, 0)
%   [P,val] = SEvalProp(Sys, P, phi, [3,7]);
%   val
%   idx_phi = find(strcmp('phi',P.props_names)); % find the index of phi in formula evaluations
%   P.props_values(idx_phi).val % fist value is equal to val
%   
%   P = SEvalProp(Sys, P, phi, 0:0.01:10);
%   PplotFormula(Sys, P, phi);
% 
%See also STL_Formula DiscrimPropValues ComputeTraj Sselect PplotFormula
%Spurge_props
%

% check arguments
if ~exist('VERBOSE','var')
    VERBOSE=0;   
end

if(~exist('method','var')||isempty(method))
    method = 'thom';
end

if ~exist('break_level','var')
    break_level = 0;
end

if(break_level>0)
    phis_tmp = [];
    for ii = 1:numel(phis)
        broken_props = STL_Break(phis(ii),break_level);
        phis_tmp = [phis_tmp broken_props(:)]; %#ok<AGROW>
    end
    phis = phis_tmp;
end

if(~exist('ipts','var')||isempty(ipts))
    ipts = 1:size(P.pts,2);
end

if(~exist('taus','var')||isempty(taus))
    taus = cell(0); % taus is not defined, we set it to empty
elseif isscalar(taus)
    tau_tmp = taus;
    taus = cell(1,numel(phis)); % create line cell array
    [taus{:}] = deal(tau_tmp); % full of taus' value
    clear('tau_tmp');
elseif isvector(taus)
    if(numel(phis)==1)
        taus = {taus};
    elseif(numel(phis)==numel(taus)) % we guess, there is one tau for each formula
        taus = reshape(taus,1,[]);
        taus = cell2mat(taus,1,ones(1,numel(taus)));
    else % all taus are for all formulas
        tau_tmp = taus;
        taus = cell(1,numel(phis));
        [taus{:}] = deal(tau_tmp);
        clear('tau_tmp');
    end
elseif(numel(taus)<numel(phis)) % it is a cell array of the wrong size
    error('SEvalProp:badTausSize','The size of taus and the size of phis are differents.');
elseif(numel(taus)<numel(phis)) % it is a cell array of the wrong size
    warning('SEvalProp:strangeTausSize','The size of taus is higher than the size of phis.');
end

if ~isfield(P,'traj')
    error('SEvalProp:noTrajField','P has no traj field.')
end
if ~isfield(P,'traj_ref')
    P.traj_ref = 1:numel(P.traj);
end
if any(P.traj_ref(ipts)==0)
    error('SEvalProp:trajNotComputed','A trajectory is not computed.');
end
if ~isfield(P,'props')
    P.props = [];
end
if ~isfield(P,'props_names')
    P.props_names = {} ;
end

val = zeros(numel(phis),numel(ipts)); %initialize array containing truth values for each property and each param set
props_values(1:size(P.pts,2)) = deal(struct()); % Temporary line containing the evaluation of the formula
for np = 1:numel(phis) % for each property
    
    phi = phis(np);  % phi = current formula
    if isa(phi, 'BreachRequirement')
        phi = STL_Formula(phi.req_monitors{1}.formula_id);
    end
    phi_name =  get_id(phi);
    i_phi = find_prop(P,phi_name);
    
    if(i_phi==0)
        % if the property does not exist in P, we add it to P
        P.props_names = [P.props_names {phi_name}];
        P.props = [P.props phi];
        i_phi = numel(P.props_names);
    end
    
    phi = STL_OptimizePredicates(Sys,phi);
    if(VERBOSE)
        fprintf('Checking  %s  on %d parameter vector(s)\n',phi_name,numel(ipts));
        fprintf('[             25%%           50%%            75%%               ]\n ');
        iprog = 0; %idx of progression bar
    end
    
    for ii = 1:numel(ipts) % we compute the truth value of phi for each parameter vector
        i_pt = ipts(ii);
        traj_tmp = P.traj(P.traj_ref(i_pt));
        Ptmp = Sselect(P,i_pt);
        
        if(isempty(taus)||isempty(taus{np}))
            [props_values(i_pt).val, props_values(i_pt).tau] = STL_Eval(Sys, phi, Ptmp, traj_tmp, traj_tmp.time(1), method);
        else
            [props_values(i_pt).val, props_values(i_pt).tau] = STL_Eval(Sys, phi, Ptmp, traj_tmp, taus{np}, method);
        end
        
        val(np,ii) = props_values(i_pt).val(1); % fill val with first evaluation at the first time point
        if isnan(val(np,ii))
            warning('SEvalProp:NaNEval','Warning: property evaluated to NaN');
        end
        
        if(VERBOSE)
            while(floor(60*ii/numel(ipts))>iprog)
                fprintf('^');
                iprog = iprog+1;
            end
        end
    end
    P.props_values(i_phi,:) = props_values; % we copy all evaluation in once to avoid inconsistent parameter set
    if(VERBOSE)
        fprintf('\n');
    end
end

end

function idx = find_prop(P,st)
%FIND_PROP finds the index of a property in a parameter set.
%
% Synopsis: idx = find_prop(P,st)
%
% Input:
%  - P  the parameter set containing the evaluation of properties
%  - st a string describing the name of the searched property
%
% Output:
%  - the index of the property evaluation if found, 0 otherwise
%

if ~isfield(P,'props_names')
    idx = 0;
    return;
else
    for idx = 1:numel(P.props_names)
        if strcmp(st,P.props_names{idx})
            return;
        end
    end
end

idx=0; % in case it is not found

end
