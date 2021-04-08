function res = tomlab_wrapper(x, Prob)

brProblem = Prob.brProblem;

res = brProblem.objective(x);

end