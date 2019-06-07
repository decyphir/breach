BrDemo.InitAutotrans;

phi = STL_Formula('phi', 'alw speed[t]<60');
B = BrAutotrans_3ranges.copy();

%% 
B.Sim();
pb = ParamSynthProblem(B, phi, 'throttle_u0', [0 100]);
pb.solve();
