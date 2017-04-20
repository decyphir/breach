%% Analysis of a simple Neural-Network based controller. 
InitBreach
narmamaglev_v1

%% Model and inputs
u_ts = 0.001;
mdl = 'narmamaglev_v1';
B = BreachSimulinkSystem(mdl); 
B.SetTime(0:.01:20);
