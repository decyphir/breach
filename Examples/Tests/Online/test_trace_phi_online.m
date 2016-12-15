function report = test_trace_phi_online(traj,phi,plot_res)
% Test online algo for one trace and one formula

if nargin ==1
    arg = traj;
    traj=arg.traj;
    phi=arg.phi;
end

if ~exist('plot_res', 'var')
    plot_res = 0;
end

Sys = CreateExternSystem('test', size(traj.X,1), 0);
P = CreateParamSet(Sys);

st_phi = disp(phi,0);
rob_phi = [inf;-inf];
ltraj  = numel(traj.time);

for it_f = 2:ltraj
%for it_f = ltraj:ltraj

%for it_f = 6:6
%    clc    
%    it_f
    traji.time = traj.time(1:it_f);
    traji.X = traj.X(:,1:it_f);
    [t, r, tl, rl, tu,ru]  =  stl_eval_mex(st_phi, [traji.time; traji.X], [0 0]);
    if (isempty(rl))
        rl= -100;
    end
    if (isempty(ru))
        ru = 100;
    end
    
    rob_phi(:,end+1)= [ru(1);  rl(1)];
    
    [rob_thom, time_thom] = STL_EvalThom(Sys, phi,P,  traji,0);
 
    status = [abs(rob_thom(1)-r(1)) < 100*eps, rl(1) <= ru(1)+100*eps, rl(1) <= r(1)+100*eps, r(1) <= ru(1)+100*eps, ...
              rob_phi(1,end)-rob_phi(1,end-1)<=100*eps, rob_phi(2,end)-rob_phi(2,end-1)>=-100*eps];
          
end

rob_max = rob_phi(1,:);
rob_min = rob_phi(2,:);

bigM = 4;
rob_max(rob_max>=100) = bigM;
rob_min(rob_min<=-100) = -bigM;

report.status= status;
report.traj = traji;
report.phi = phi;
report.thom.time = time_thom;
report.thom.values = rob_thom;
report.online.rob.time = t;
report.online.rob.values = r;
report.online.rob_low.time = tl;
report.online.rob_low.values = rl;
report.online.rob_up.time = tu;
report.online.rob_up.values = ru;
report.online.rob_phi = rob_phi;
% status = [rob_thom == rob, rob_low <= rob_up, rob_low <= rob, rob <= rob_up]

if plot_res
    figure
    hold on;
    grid on;
    eps_ = .01;
    
    plot(traj.time, traj.X(1,:),'b',traj.time, traj.X(2,:),'g', traj.time, 0*traj.X(1,:), 'y');
    plot(traj.time, 0*traj.X(1,:)+rob_thom(1), 'r--');
    plot(traj.time, 0*traj.X(1,:)+r(1), 'g--');
    legend({'x1', 'x2', 'zero', 'rob thom', 'rob'});
    %plot(traj.time, traj.X(1,:),'b', traj.time, 0*traj.X(1,:)-bigM-eps_, 'k--');
        
    text(0.1, bigM+eps_-.1, '+\infty', 'FontSize', 24);
    text(0.1, -bigM-eps_+.2, '-\infty','FontSize', 24);
    
    stairs(traji.time(1:end),rob_max,'r');
    stairs(traji.time(1:end),rob_min, 'k')
    title(st_phi)
end


end

