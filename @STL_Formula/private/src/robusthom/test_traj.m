function traj = test_traj(n)

traj1.time = 1:n;
traj1.X = 4*rand([1 n])-2;

traj2.time = 1:n;
traj2.X = 4*rand([1 n])-2;

traj.time = linspace(0, 100 ,n);
traj.X = [ traj1.X(1,:); traj2.X(1,:) ];
traj.param = [0 0 0 0 0 0];

