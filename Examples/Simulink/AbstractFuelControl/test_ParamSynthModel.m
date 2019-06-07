BrDemo.InitAutotrans;
phi = STL_Formula('phi', 'alw speed[t]<60');
pb = ParamSynthProblem(B, phi, 'throttle_u0', [0 100]);
pb.display = 'on';
pb.setup_global_nelder_mead();
pb.solve();
