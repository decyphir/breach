function [val, vald] = generic_predicate(~, traj, t, params) 
%GENERIC_PREDICATE implements a generic predicate
% 
% Synopsis: [val, vald] = generic_predicate(mode, traj, t, params)
% 
% Inputs:
%  - mode   : some switch variable (not used) 
%  - traj   : one trajectory. Must contain fields time and X. May
%             contains field XS.
%  - t      : time values in which we evaluate the predicate 
%  - params : structure which must have fields Sys (the system), P (the
%             parameter set) and fn (the evaluation function). It may
%             contain fields pre_pred (which contains fields init_par_st,
%             init_var_st and eval_fn) and eval_pred.
% 
% Outputs:
%  - val  : values of the predicate at times. NaN if there is a fatal
%           error while evaluating the formula.
%  - vald : derivatives of the values (not used)
% 
%See also STL_Formula
%

if isfield(params,'eval_pred')
    params.eval_pred(params.P);
end
    
global BreachGlobOpt;
eval(BreachGlobOpt.GlobVarsDeclare);

if numel(traj.time)==1
    dt__ = 1;    %  not sure this makes sense, will see
else
    dt__ = traj.time(1,2)-traj.time(1,1);
end
Sys = params.Sys;
P =params.P;

if isfield(params, 'pre_pred') % if the predicate has been optimized
    pre_pred = params.pre_pred;
       
    eval(pre_pred.init_par_st);
    eval(pre_pred.init_var_st);
    
    val = zeros(1, numel(t));
    val(1,:) = eval(pre_pred.eval_fn);
    if isfield(BreachGlobOpt, 'RobustSemantics') 
      switch (BreachGlobOpt.RobustSemantics)
          case -1
              val = ltr(val,dt__);
          case 1
              val = rtr(val, dt__);
          case 2
              rval = abs(rtr(val,dt__));
              lval = abs(ltr(val, dt__));
              sval = sign(val);
              val = sval.*min([rval; lval]);
      end
      if isfield(BreachGlobOpt, 'NormalizePredicates')
          if BreachGlobOpt.NormalizePredicates
              val = normalize_rob(val);
          end
      end
    end
    return;
end

fn_ = params.fn;

% don't like that... But not so trivial to remove
traj.time = [-1e99 traj.time 1e99];          
traj.X = [traj.X(:,1) traj.X traj.X(:,end)];  
 
if isfield(traj,'XS')
    traj.XS = [traj.XS(:,1) traj.XS traj.XS(:,end)];
end

% we need the values of the parameters in the workspace

fn_ = regexprep(fn_,'([^\.])/', '$1\./');
fn_ = regexprep(fn_,'([^\.])\*', '$1\.*');
fn_ = regexprep(fn_,'([^\.])\^', '$1\.^');

% we capture variables involved in the formula
[~,~,~,~,tokens] = regexp(fn_, ['^(\w+)\[.+?\]|'  ...
    '\W(\w+)\[.+?\]|' ...
    'ddt\{\s*(\w+)\s*\}\[.+?\]|' ...
    'd\{\s*(\w+)\s*\}{.+?}\[.+?\]'...
    '\{(.+?)\}[.+?]']);

SigList = unique(cat(2,tokens{:}));
for ii_var = 1:numel(SigList)
    pcurr = SigList{ii_var};
    i_var = FindParam(Sys, pcurr);
        
    % test if current variable is found at the beginning of the predicate
    [start_idx, end_idx, ~, matches, tokens] = regexp(fn_, ['^' pcurr '\[(.+?)\]']);
    
    if numel(matches)
        pcurr = tokens{1}{1};
        time = eval(tokens{1}{:});
        newpcurr = sprintf('%s__LB__0_RB_',pcurr); %NM: Nota: index 1 may be used later, so we use 0 !
        X = interp1(traj.time, traj.X(i_var,:), time, 'linear', traj.X(i_var,end)); %#ok<NASGU>
        eval(sprintf('%s = X;',newpcurr)); % faster than call interp1 in eval
        fn_ = [fn_(1:start_idx-1) newpcurr fn_(end_idx+1:end)];
    end
    
    
    % current variable inside curly brackets
    [start_idx, end_idx, ~, matches, tokens] = regexp(fn_, ['\{' pcurr '\}' '\[(.+?)\]']);
    
    if numel(matches)
        pcurr = regexprep(tokens{1}{1}, '\W', '__');
        time = eval(tokens{1}{:});
        newpcurr = sprintf('%s__LB__0_RB_',pcurr); %NM: Nota: index 1 may be used later, so we use 0 !
        X = interp1(traj.time, traj.X(i_var,:), time, 'linear', traj.X(i_var,end)); %#ok<NASGU>
        eval(sprintf('%s = X;',newpcurr)); % faster than call interp1 in eval
        fn_ = [fn_(1:start_idx-2) newpcurr fn_(end_idx+2:end)];
    end
    
    
    % test if current variable is found inside the predicate
    [start_idx, end_idx, ~, matches, tokens] = regexp(fn_, ['\W' pcurr '\[(.+?)\]']); % \W makes sure we find its complete name (not just suffix)
    
    start_idx = start_idx+1; % skips \W
    for jj = 1:numel(matches)
        time = eval(tokens{jj}{:});
        newpcurr = sprintf('%s__LB__%d_RB_',pcurr,jj);
        X = interp1(traj.time, traj.X(i_var,:), time, 'linear', traj.X(i_var,end)); %#ok<NASGU>
        eval(sprintf('%s = X;',newpcurr)); % faster than call interp1 in eval
        fn_ = [fn_(1:start_idx(jj)-1) newpcurr fn_(end_idx(jj)+1:end)]; % replace xx[tt] by xx__LB__jj_RB_
        lj = end_idx(jj)-start_idx(jj);
        start_idx = start_idx-lj-1+length(newpcurr); % update start_idx and end_idx so they
        end_idx = end_idx-lj-1+length(newpcurr);     % still correspond to begin and end of matches
    end
    
    % time derivative
    [start_idx, end_idx, ~, matches, tokens] = regexp(fn_, ['ddt\{\s*' pcurr '\s*\}\[(.+?)\]']);
    
    for jj = 1:numel(matches)
        time = eval(tokens{jj}{:});
        newpcurr = sprintf('%s__DLB__%d_DRB_',pcurr,jj);
        X = interp1(traj.time, traj.X', time , 'linear', NaN)'; % X is a column of size Sys.DimX x 1 ; NaN = temp value
        X(isnan(X)) = traj.X(isnan(X),end); % replace NaN by the last value of the traj
        Pts = traj.param';    % Pts has size Sys.DimP x 1
        Pts(1:Sys.DimX,:) = X;
        Fval = GetF(Sys, time, Pts); % ask for derivatives at time instant "time"
        X = Fval(i_var,:); %#ok<NASGU>
        eval(sprintf('%s = X;',newpcurr));
        fn_ = [fn_(1:start_idx(jj)-1) newpcurr fn_(end_idx(jj)+1:end)];
        lj = end_idx(jj)-start_idx(jj);
        start_idx = start_idx-lj-1+length(newpcurr); % update start_idx and end_idx so they
        end_idx = end_idx-lj-1+length(newpcurr);  % still correspond to begin and end of matches
    end
    
    % parameter derivative
    [start_idx, end_idx, ~, matches, tokens] = regexp(fn_, ['d\{\s*' pcurr '\s*\}{(.+?)}\[(.+?)\]']);
    
    for jj = 1:numel(matches)
        time = eval(tokens{jj}{2});
        sparam = tokens{jj}{1};
        newpcurr = sprintf('Sensi_%s__DLB__%s_DRB_',pcurr,sparam);
        
        i_param = FindParam(params.P, sparam);
        [xs, traj] = GetTrajSensi(Sys, traj, i_var, i_param);
        X = interp1(traj.time, xs, time, 'linear', xs(end)); %#ok<NASGU>
        eval(sprintf('%s = X;',newpcurr));
        fn_ = [fn_(1:start_idx(jj)-1) newpcurr fn_(end_idx(jj)+1:end)];
        lj = end_idx(jj)-start_idx(jj);
        start_idx = start_idx-lj-1+length(newpcurr);
        end_idx = end_idx-lj-1+length(newpcurr);
    end
end

val = eval(fn_);

if isscalar(val)
    val = val*ones(1, numel(t));
end

if isfield(BreachGlobOpt, 'RobustSemantics')
    switch (BreachGlobOpt.RobustSemantics)
        case -1
            val = ltr(val, dt__);
        case 1
            val = rtr(val, dt__);
        case 2
            rval = abs(rtr(val, dt__));
            lval = abs(ltr(val,dt__));
            sval = sign(val);
            val = sval.*(rval+lval)/2;                
    end            
end
if isfield(BreachGlobOpt, 'NormalizePredicates')
    val = normalize_rob(val);        
end

    
end


function val = normalize_rob(val)
    % normalize values into the range [-1, +1] by scaling by 
    % max(abs(min),abs(max))
    scaling_factor = max(abs(min(val)),  abs(max(val)));
    if scaling_factor == 0
        % If all values are 0, scaling_factor = 0
        % Avoid division by zero by just returning val as it is
        return
    end
    val = val/scaling_factor;   
end