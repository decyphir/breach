% mex RobustEv.cpp robustness.cpp signal.cpp mex_routines.cpp
traj = test_traj(20);

t = traj.time;
y1 = traj.X(1,:);
y2 = traj.X(2,:);

[to yevinf] = RobustEv(t,y1, [0 inf]);
figure
plot( t, y1, to,yevinf )
legend({'y','ev y'})

[to yevainf] = RobustEv(t,y1, [5. inf]);
figure
plot( t, y1, to,yevainf )
legend({'y','ev y'})

[to yevb] = RobustEv(t,y1, [0 12.]);
figure
plot( t, y1, to,yevb )
legend({'y','ev y'})

[to yevab] = RobustEv(t,y1, [10. 20.]);
figure
plot( t, y1, to, yevab )
legend({'y','ev y'})
