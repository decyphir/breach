Init_BrAFC;

BrAFC.Sim();
BrAFC.UpdateSignalRanges();

phi_in_range = STL_Formula('phi_in_range','alw (DistFromRange(traj,P.SigRange,t)>-0.01)');
phi_in_range= set_params(phi_in_range, 'SigRange',BrAFC.SignalRanges);

params = {'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}
param_ranges = [10 20; 10 60]

prob = FalsificationProblem(BrAFC, phi_in_range, params, param_ranges );

prob.solver_options.nb_samples = 50;
prob.solve();
fval = prob.obj_best;


while fval<0
    
    % get and computes the falsifying trace
    BrFalse = prob.GetBrSet_False();
    BrFalse.Sim();
    
    % Concat the trace to previous ones and update ranges
    BrAFC.Concat(BrFalse)
    BrAFC.UpdateSignalRanges();
    
    % Update formula to falsify
    phi_in_range= set_params(phi_in_range, 'SigRange',BrAFC.SignalRanges);
    
    % Update Falsification problem
    prob = FalsificationProblem(BrFalse, phi_in_range ,  params, param_ranges);
    prob.solver_options.nb_samples = 50;
    prob.solve();
    fval = prob.obj_best;

end