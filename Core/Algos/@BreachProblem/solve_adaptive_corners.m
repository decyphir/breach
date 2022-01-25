function res = solve_adaptive_corners(this, varargin)

BrC = BreachSet(this.params);
for ip = 1:numel(this.params)
    BrC.SetDomain(this.params{ip}, this.BrSys.GetDomain(this.params{ip}));
end

num_corners =  min(this.max_obj_eval,this.solver_options.num_corners);

% Creates BreachSet with one param per group
sigs = this.BrSys.InputGenerator.GetAllSignalsList();
grouparams = this.params;
groupinputs_reps = {};
ig=1;
for is = 1:numel(sigs)
    groupin = intersect(this.BrSys.expand_param_name([sigs{is} '_']), this.params);
    if ~isempty(groupin)
        groupinputs{ig} = groupin;
        groupinputs_reps{ig} = [sigs{is} '_group'];
        grouparams = setdiff(grouparams, groupinputs{ig});
        ig = ig+1;
    end
end
if isempty(grouparams)
    BrC_group = BreachSet(groupinputs_reps);
else
    BrC_group = BreachSet([grouparams groupinputs_reps]); % adds remaining params to group params
end

% Sample group BreachSet
gparams= BrC_group.GetParamList();
for ig = 1:numel(gparams)
    BrC_group.SetDomain(gparams{ig}, BreachDomain('enum', [1 2], [1 2]));
end
BrC_group.CornerSample(num_corners);

% Extract concrete corners
X0 = zeros(numel(this.params), size(BrC_group.P.pts, 2));
if ~isempty(grouparams)  % TEST ME!
    for ig = 1:numel(grouparams)
        idx = strcmp(this.params, grouparams{ig});
        which_corn = BrC_group.GetParam(grouparams{ig});
        corn_vals =  [this.lb(idx) this.ub(idx)];
        X0(idx,:) =  arrayfun(@(c)(corn_vals(c)), which_corn);
    end
end
for ig = 1:numel(groupinputs_reps)
    which_corn = BrC_group.GetParam(groupinputs_reps{ig});
    for ip = 1:numel(groupinputs{ig})
        idx = strcmp(this.params, groupinputs{ig}{ip});
        corn_vals =  [this.lb(idx) this.ub(idx)];
        X0(idx,:) =  arrayfun(@(c)(corn_vals(c)), which_corn);
    end
end
if size(X0,2)<num_corners % complete with normal corners
    BrC.CornerSample(num_corners);
    X0 = unique([X0 BrC.GetParam(this.params)]', 'rows','stable')';
    X0 = X0(:, 1:min(num_corners, size(X0,2)));
end

if isfield(this.BrSys.Sys, 'use_parallel') && this.BrSys.Sys.use_parallel
    error('TODO: Parallel computation');
    
else
    % Go through each sample one at a time
    samplesSinceLastFalsification = 0;
    idxFalsified = []; % Keep track of which reqs have been falsified
    allFval = [];
    allCval = [];
    for xCounter = 1:num_corners
        [fval, cval] = this.objective(X0(:, xCounter));
        
        % Add fval, cval to their structures
        allFval(:, end+1) = fval;
        allCval(:, end+1) = cval;
        
        newIdxFalsified = find(fval < 0);
        differentFalsified = setdiff(newIdxFalsified, idxFalsified);
        if isempty(differentFalsified)
            % No new falsified reqs
            samplesSinceLastFalsification = samplesSinceLastFalsification + 1;
        else
            % New falsified reqs
            samplesSinceLastFalsification = 0;
            idxFalsified = unique([idxFalsified newIdxFalsified']);
            
            fprintf(['Iter ' num2str(xCounter) ': ' ...
                num2str(numel(differentFalsified)) ' new falsified reqs\n']);
        end
        
        % Break if we went too long without a new falsification
        if samplesSinceLastFalsification > this.solver_options.relative_threshold*num_corners
            fprintf(['Iter ' num2str(xCounter) ': ' ...
                'No new falsified reqs for ' ...
                num2str(samplesSinceLastFalsification) ...
                ' samples - stopping adaptive_corners falsification ...\n']);
            break
        end
        
    end
    
    % Create a return structure similar to solve_corners.m
    [fbest, ibest] = min(min(allFval));
    res = struct('BrSys', this.BrSys, 'X0',X0(:, 1:xCounter),'x',X0(:,ibest),'f', fbest, 'fval', allFval,'cval',allCval);
    
end


end