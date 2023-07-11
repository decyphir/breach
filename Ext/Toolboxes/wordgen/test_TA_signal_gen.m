%% TA based signal generator setup 
%
init_TA_signal_gen;

S1 = S0.copy();
S1.SetParamRanges([pevts pbranching], [0 1]);
pvar = S1.GetVariables();
S1.SampleDomain(pvar,100);

S1.Sim();
figure;
S1.PlotSignals();

%% Checking reachable labels
STL_ReadFile('Autotrans_spec.stl');
Rreach_labels = BreachRequirement({req_a, req_b, req_c, req_d, req_e, req_f, req_g, req_h});    
Rreach_labels.Eval(phi1);
BreachSamplesPlot(Rreach_labels);

%% 
S2 = S0.copy();
mdl = 'Autotrans_wordgen';
Ba = BreachSimulinkSystem(mdl);
Ba.SetTime(time);
Ba.SetInputGen(S2);
Ba.SetParamRanges([pevts pbranching], [0 1]);
Ba.SampleDomain(pvar,10);
Ba.Sim();


Rreach_labels = BreachRequirement({req_a, req_b, req_c, req_d, req_e, req_f, req_g, req_h});    
Rreach_labels.Eval(Ba)

%%
BreachSamplesPlot(Rreach_labels);

