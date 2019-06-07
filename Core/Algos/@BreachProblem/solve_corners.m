function res = solve_corners(this)
% solve_corners works with quasi-random sampling   
% 


BrC = this.BrSet.copy();
BrC.ResetParamSet();
BrC.SetParamRanges(this.params, [this.lb this.ub])
num_corners =  this.solver_options.num_corners;
BrC.CornerSample(num_corners);
X0 = BrC.GetParam(this.params);

if ~strcmp(this.display,'off')
    fprintf('\n Running %g corners\n', size(X0, 2));
    this.display_status_header();
end

res = this.FevalInit(X0);
this.add_res(res); 

end
