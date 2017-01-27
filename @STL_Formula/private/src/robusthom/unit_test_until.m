mex RobustUntil.cpp robustness.cpp signal.cpp mex_routines.cpp
traj = test_traj(20);

t = traj.time;
y1 = traj.X(1,:);
y2 = traj.X(2,:);

[to yuntinf] = RobustUntil(t,y1, t, y2, [0 inf]);
figure
plot( t, y1, t,y2, to,yuntinf )
legend({'y1', 'y2',' y1 until y2'})
%return

[to yuntainf] = RobustUntil(t,y1, t, y2, [5. inf]);
figure
plot( t, y1, t, y2, to,yuntainf )
legend({'y1', 'y2',' y1 until[5 inf] y2'})

[to yuntb] = RobustUntil(t,y1, t, y2, [0 12.]);
figure
plot( t, y1, t,y2, to,yuntb )
legend({'y1', 'y2',' y1 until[0 12] y2'})


[to yuntab] = RobustUntil(t,y1, t,y2, [10. 20.]);
figure
plot( t, y1, t, y2, to, yuntab )
legend({'y1', 'y2',' y1 until[10,20] y2'})
