%% Init Randon number generator
seed=0;
s = RandStream('mcg16807','Seed',seed);
RandStream.setGlobalStream(s);

%% generates formulas
nphis = 20; % number of formulas
sz_phis= 3; % size of formulas in number of operators

% Init cell array of fomulas 
PHIS = cell(1,nphis);
for iphi=1:nphis
    st_phi =  rand_formula(sz_phis, {'evI', 'alwI', 'and', 'or'} );
    PHIS{iphi} = STL_Formula(['phi_rand' num2str(iphi)], st_phi);
end

%% generate trajectories with variable number of samples
n_dim = 2*sz_phis;
n_traces = 10; 
n_samples = 100;

B = BreachTraceSystem(n_dim);
B.AddRandomTraces(n_traces, n_samples); 

%% test all and record results
res = []; 
fprintf('\n');
for iphi=1:nphis
    fprintf('.');
    phi = PHIS{iphi};
    B.AddSpec(phi);
    B.CheckSpec(phi);
end
fprintf('\n');