function traj = gen_traj(n,k)

traj.time = linspace(0,100,n);
traj.X = 4*rand([k n])-4*rand();
traj.param = zeros(1,k+1);

end
