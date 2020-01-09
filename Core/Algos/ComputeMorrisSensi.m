function res = ComputeMorrisSensi(R, B, varargin)
%% ComputeMorrisSensi compute sensitivities using Morris method for all requirement in R. Assumes B is a set with ranges. 
%  
%  Note: if B is a BreachSimulinkSystem with no traces, we compute the
%  Morris samples and corresponding traces and stores them in B. If B
%  already contains Morris samples and traces they are not re-computed.

opt = struct('num_path', 100, ...      % number of paths, i.e., set of
                    ...                % samples providing 1 pair of samples per dim
                'rand_seed', 1, ...    % random seed for reproducibility    
                'size_grid', 5 ...     % number of grid levels, intervals in each dim
    );

opt = varargin2struct(opt, varargin{:});

[vars, iv] = B.GetVariables(); % variables and indices
    
if ~isfield(B.P, 'opt_morris')||~isequal(opt, B.P.opt_morris)  
    ranges = B.GetParamRanges(iv);    
    Sys= CreateSystem({},vars, ranges(:,1));
    P  = CreateParamSet(Sys, vars, ranges);        
    r = opt.num_path;
    p = opt.size_grid;
    s = opt.rand_seed;
    Pr = pRefine(P, p,r,s);
    X0 = Pr.pts;
    B.ResetParamSet;
    B.SetParam(vars, X0)
    B.P.opt_morris =opt;
    B.P.D = Pr.D;       
end

R.Eval(B);

for ir =  1:numel(R.req_monitors)
    res{ir}.params = vars;
    res{ir}.rob = R.traces_vals(:,ir)';
    [res{ir}.mu, res{ir}.mustar, res{ir}.sigma, res{ir}.sigmastar, res{ir}.EE] = EEffects(res{ir}.rob, B.P.D, opt.size_grid);
    fprintf(['\nSensitivities for ' R.req_monitors{ir}.name]); 
    display_morris_result(res{ir});
end

end