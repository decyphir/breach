

%% Init Breach and compile mex functions
clear
InitBreach;
compile_stl_mex;

%% Always use the same random numbers
seed=0;
s = RandStream('mcg16807','Seed',seed);
RandStream.setGlobalStream(s);

%% generates formulas
nphis = 1;
sz_phis= 1;

% defines predicates and formula
for i = 1:2*sz_phis
    STL_Formula( ['mu' num2str(i)] , ['x' num2str(i) '[t]>0']);
end

PHIS = cell(1,nphis);
for iphi=1:nphis
    st_phi =  rand_formula(sz_phis, {'evI', 'alwI', 'and', 'or'} );
    PHIS{iphi} = STL_Formula('phi', st_phi);
end

%% generate trajectories with variable number of samples
n_samples = [20 100];
for it = 1:numel(n_samples)
    trajs(it) = gen_traj(n_samples(it), 2*sz_phis);
end

report_bug = [];
for iphi=1:nphis
    phi = PHIS{iphi}
    for it = 1:numel(trajs)
        reports(it,iphi)= test_trace_phi_online(trajs(it), phi);%,1); drawnow
        if ~all(reports(it,iphi).status(1:end))
            if isempty(report_bug)
                report_bug = reports(it,iphi);
            else
                report_bug(end+1)= reports(it,iphi);               
            end
            warning(interpret_status(report_bug(end)));
        end
    end
end

