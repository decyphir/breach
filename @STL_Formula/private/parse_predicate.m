function mu = parse_predicate(S, mu)
%PARSE_PREDICATE pre_parses a predicate to make its evaluation faster
% 
% Synopsis: pre_pred = parse_predicate(S, mu)
% 
% Inputs:
%   - S  : can be a system Sys or a parameter set P.
%   - mu : the predicate to optimize
% 
% Output:
%  - mu : the optimized predicate. This method adds a field pre_pred to the
%         field params of mu. This field pre_pred contains three fields
%         init_var_st, init_par_st (which contain string to execute before
%         the evaluation of the predicate) and eval_fn. They then avoid to
%         parse the predicate each time it is evaluated.
% 
%See also STL_OptimizePredicates
%

fn_ = mu.params.fn;

% replaces * and
fn_ = regexprep(fn_,'([^\.])/', '$1\./');
fn_ = regexprep(fn_,'([^\.])\*', '$1\.*');
fn_ = regexprep(fn_,'([^\.])\^', '$1\.^');

pre_pred.init_par_st = '';
pre_pred.init_var_st = '';

for ii=1:S.DimX
    pcurr = S.ParamList{ii};
    
    % test if current variable is found at the beginning of the predicate
    [start_idx, end_idx, ~, matches, tokens] = regexp(fn_, ['^' pcurr '\[(.+?)\]']);
          
    for jj = 1:numel(matches)
        time = tokens{jj}{:}; % time is then the string indicating at which time the variable is evaluated
        newpcurr = [pcurr '__LB__' int2str(jj) '_RB_'];
        % CAUTION: this could lead to an error as __LB_1_RB_ may be defined
        % twice. This does not happend thanks to the way the evaluation
        % function is defined (no function starts with a predicate name)
        pre_pred.init_var_st = [pre_pred.init_var_st [newpcurr '= interp1(traj.time, traj.X(' int2str(ii) ',:),' time ',''linear'',' 'traj.X(' int2str(ii) ',end)'  ');' char(10)]];
        fn_ = [fn_(1:start_idx(jj)-1) newpcurr fn_(end_idx(jj)+1:end)];
        lj = end_idx(jj)-start_idx(jj);
        start_idx = start_idx-lj-1+length(newpcurr);
        end_idx = end_idx-lj-1+length(newpcurr);
    end
    
    % test if current variable is found inside the predicate
    [start_idx, end_idx, ~, matches, tokens] = regexp(fn_, ['\W' pcurr '\[(.+?)\]']); % \W makes sure we find its complete name (not just suffix)
    
    start_idx = start_idx+1;
    for jj = 1:numel(matches)
        time = tokens{jj}{:};
        newpcurr = [pcurr '__LB__' int2str(jj) '_RB_'];
        pre_pred.init_var_st = [pre_pred.init_var_st newpcurr '= interp1(traj.time, traj.X(' int2str(ii) ',:),' time ',''linear'',' 'traj.X(' int2str(ii) ',end)' ');' char(10)];
        %      eval ([newpcurr '= interp1(traj.time, traj.X(' num2str(i) ',:), time,''linear'',inf);']);
        fn_ = [fn_(1:start_idx(jj)-1) newpcurr fn_(end_idx(jj)+1:end)];
        lj = end_idx(jj)-start_idx(jj);
        start_idx = start_idx-lj-1+length(newpcurr);
        end_idx = end_idx-lj-1+length(newpcurr);
    end
    
    [start_idx, end_idx, ~, matches, tokens] = regexp(fn_, ['ddt\{\s*' pcurr '\s*\}\[(.+?)\]']);
    
    for jj = 1:numel(matches)
        time = tokens{jj}{:};
        newpcurr = [pcurr '__DLB__' int2str(jj) '_DRB_'];
        pre_pred.init_var_st = [pre_pred.init_var_st ['X = interp1(traj.time, traj.X'','  time ',''linear'',NaN)'';' char(10)]];
        pre_pred.init_var_st = [pre_pred.init_var_st ['X(isnan(X)) = traj.X(isnan(X),end);' char(10)]];
        pre_pred.init_var_st = [pre_pred.init_var_st ['Pts = repmat(traj.param'', [1 size(X,2)]);' char(10) ]];
        pre_pred.init_var_st = [pre_pred.init_var_st ['Pts(1:Sys.DimX,:) = X;'  char(10)]];
        pre_pred.init_var_st = [pre_pred.init_var_st ['Fval = GetF(Sys,' time ', Pts);'  char(10)]];
        pre_pred.init_var_st = [pre_pred.init_var_st [newpcurr ' = Fval(' int2str(ii) ',:);'] char(10)];
        
        fn_ = [fn_(1:start_idx(jj)-1) newpcurr fn_(end_idx(jj)+1:end)];
        lj = end_idx(jj)-start_idx(jj);
        start_idx = start_idx-lj-1+length(newpcurr);
        end_idx = end_idx-lj-1+length(newpcurr);
    end
    
    [start_idx, end_idx, ~, matches, tokens] = regexp(fn_, ['d\{\s*' pcurr '\s*\}{(.+?)}\[(.+?)\]']);
    
    for jj = 1:numel(matches)
        time = tokens{jj}{2};
        sparam = tokens{jj}{1};
        
        is_ = FindParam(S,sparam);
        newpcurr = ['Sensi_' pcurr '__DLB__' sparam '_DRB_'];   
        pre_pred.init_var_st = [pre_pred.init_var_st ['[xs_, traj] = GetTrajSensi(Sys, traj,' int2str(ii) ',' int2str(is_) ' );' char(10)]];
        %NM : Pay attention here : does it really work?! (use of an array as extravalp ; see above for a potential fix)
        pre_pred.init_var_st = [pre_pred.init_var_st [newpcurr '= interp1(traj.time, xs_'',' time ',''linear'',xs_(:,end)'');' char(10)]];
        
        fn_ = [fn_(1:start_idx(jj)-1) newpcurr fn_(end_idx(jj)+1:end)];
        lj = end_idx(jj)-start_idx(jj);
        start_idx = start_idx-lj-1+length(newpcurr);
        end_idx = end_idx-lj-1+length(newpcurr);
    end
end

pre_pred.eval_fn = fn_;
mu.params.pre_pred = pre_pred;

end

