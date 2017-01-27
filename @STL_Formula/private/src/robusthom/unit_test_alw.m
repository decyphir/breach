% mex RobustEv.cpp robustness.cpp signal.cpp mex_routines.cpp
traj = test_traj(10000000);

t = traj.time;
y1 = traj.X(1,:);
y2 = traj.X(2,:);

tp = [-1e100 t 1e100];
yp = [y1(1) y1 y1(end)];

i2 = 20;
t= [t t(end)+i2];
y1 = [y1 y1(end)];
[to yalwinf] = RobustEv(t,-y1, [0. i2]);

return;
figure
plot( t, y1, to,-yalwinf )
legend({'y','alw y'})

% [to yalwainf] = RobustEv(t,-y1, [5. inf]);
% figure
% plot( t, y1, to,-yalwainf )
% legend({'y','ev y'})
% 
% [to yalwb] = RobustEv(t,-y1, [0 12.]);
% figure
% plot( t, y1, to,-yalwb )
% legend({'y','alw y'})
% 
% [to yalwab] = RobustEv(t,-y1, [10. 20.]);
% figure
% plot( t, y1, to, -yalwab )
% legend({'y','alw y'})
