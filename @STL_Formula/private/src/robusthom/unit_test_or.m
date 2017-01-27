traj = test_traj(10000000);

t = traj.time;
y1 = traj.X(1,:);
y2 = traj.X(2,:);

[to yo] = RobustOr(t,y1,t,y2);
return;
figure
plot( t, y1, t, y2 , to,yo );
legend({'y1', 'y2','y1 and y2'});
