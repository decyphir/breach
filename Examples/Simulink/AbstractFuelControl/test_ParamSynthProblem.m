Init_BrAFC;
AFCp = BrAFC.copy();
AFCp.SetParamRanges({'Pedal_Angle_pulse_period','Pedal_Angle_pulse_amp'},  [10 20; 5 60]);
AFCp.QuasiRandomSample(10);
AFCp.Sim();

STL_ReadFile('simple_spec.stl');

% will find tol for all trajectories in AFCp
synth_problem = ParamSynthProblem(AFCp, AF_alw_ok, {'tol'}, [0 0.1]);
synth_problem.setup_solver('fmincon')
synth_problem.solve();

% Let's check that the best tol found works:
tol_best = synth_problem.x_best
rob = synth_problem.robust_fn(tol_best) % get robustness for each 10 trajectories
                                        % -> all positive, one tight (close
                                        % to 0)