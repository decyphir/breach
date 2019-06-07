clear;
InitBreach;
 
M = csvread('Data_2.csv');
a = M(:,[1:2]);
b = M(:,[1,3]);
t = transpose(M(:,1));
 
% Interface workspace data with a Breach object
workspace_data = from_workspace_signal_gen({'a','b'});
B = BreachSignalGen({workspace_data});
B.Sim(t);

%% EXAMPLE 1 -----------------------

% define the formula
phi1 = STL_Formula('phi1', 'alw(a[t] <= 3)');
phi1 = set_out_signal_names(phi1, {'a', 'b'});

% Set up the Breach requirement
MyReq1 = BreachRequirement(phi1);
MyReq1.Eval(B);
verdict = MyReq1.Explain(B, phi1);
MyReq1.PlotDiag_debug(verdict);

%% EXAMPLE 2 -----------------------

% define the formula
phi2 = STL_Formula('phi2', 'alw(a[t] <= 5)');
phi2 = set_out_signal_names(phi2, {'a', 'b'});

% Set up the Breach requirement
MyReq2 = BreachRequirement(phi2);
verdict = MyReq2.Explain(B, phi2);
MyReq2.PlotDiag(phi2, verdict);

%% EXAMPLE 3 -----------------------

% define the formula
phi3 = STL_Formula('phi3', 'alw((a[t] <= 4 and a[t] > 3) => ev_[1,3] (b[t] <= 3))');
phi3 = set_out_signal_names(phi3, {'a', 'b'});

% Set up the Breach requirement
MyReq3 = BreachRequirement(phi3);
verdict = MyReq3.Explain(B, phi3);
MyReq3.PlotDiag(phi3, verdict);

%% EXAMPLE 4 -----------------------

% define the formula
phi4 = STL_Formula('phi4', 'alw((a[t] >= 4) => ev_[1,2] (b[t] >= 6))');
phi4 = set_out_signal_names(phi4, {'a', 'b'});

% Set up the Breach requirement
MyReq4 = BreachRequirement(phi4);
verdict = MyReq4.Explain(B, phi4);
MyReq4.PlotDiag(phi4, verdict);

%% EXAMPLE 5 -----------------------

% define the formula
phi5 = STL_Formula('phi5', 'alw((a[t] >= 2) => ev_[1,2] (b[t] <= -1))');
phi5 = set_out_signal_names(phi5, {'a', 'b'});

% Set up the Breach requirement
MyReq5 = BreachRequirement(phi5);
verdict = MyReq5.Explain(B, phi5);
MyReq5.PlotDiag(phi5, verdict);

%% EXAMPLE 6 -----------------------

% define the formula
phi6 = STL_Formula('phi6', 'alw((a[t] >= 2) => ev_[1,2] (b[t] <= 10))');
phi6 = set_out_signal_names(phi6, {'a', 'b'});

% Set up the Breach requirement
MyReq6 = BreachRequirement(phi6);
verdict = MyReq6.Explain(B, phi6);
MyReq6.PlotDiag_debug(verdict);

%% EXAMPLE 7 -----------------------

% define the formula
phi7 = STL_Formula('phi7', 'alw((not (a[t] == 0)) => (b[t] == 0))');
phi7 = set_out_signal_names(phi7, {'a', 'b'});

% Set up the Breach requirement
MyReq7 = BreachRequirement(phi7);
verdict = MyReq7.Explain(B, phi7);
MyReq7.PlotDiag_debug(verdict);


