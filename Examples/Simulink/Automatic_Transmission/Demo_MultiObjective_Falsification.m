BrDemo.InitAutotrans

%% Multi-requirement 
% 
req1 = STL_Formula('phi_speed', 'not (alw_[35,40] (speed[t]<speed_up and speed[t]>speed_low))'); 
req2 = STL_Formula('phi_rpm', 'not (alw_[35,40] (RPM[t]<rpm_up and RPM[t]>rpm_low))');

R = BreachRequirement({req1, req2});
R.SetParam({'speed_low','speed_up'}, [80, 85]);
R.SetParam({'rpm_low','rpm_up'}, [3000, 3100]);

pb = FalsificationProblem(BrAutotrans_3ranges, R);
pb.SetupParallel();
pb.solve();

%%
%
Rlog = pb.GetLog();
BreachSamplesPlot(Rlog);

