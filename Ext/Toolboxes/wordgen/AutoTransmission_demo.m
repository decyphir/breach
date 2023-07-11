%%% Old Script for demo of Classification based method 
bdclose all
close all
clear
InitBreach


%% Breach Interface Object Creation

model_name = 'autotrans_mod04';
fprintf('\n Creating breach interface with simulink model %s\n',model_name)

simTime = 30 ; 
fprintf('\n Simulation time horizon is %d seconds\n',simTime)

fprintf('\n Press any key to continue')
pause

BrSys = CoverageBreachSet(model_name,{});
BrSys.SetTime([0 simTime]);

%% Set input signals

fprintf('\n Parametrizing input signals throttle and break....\n')
fprintf('\n Input signals are parametrized as piecewise constant\n ')
input_gen.type = 'UniStep';

N1 = 7; N2 = 3;
fprintf('Number of control points for throttle input is %d\n',N1)
fprintf('Number of control points for break input is %d\n',N2 )
input_gen.cp = [N1 N2];
input_gen.method = {'previous','previous'};
BrSys.SetInputGen(input_gen);

fprintf('\n Press any key to continue\n')
pause

%% Specifying parameter names
for i=0:N1-1
    signal_u0{1,i+1}=strcat('In1_u',num2str(i));
end

for i=0:N2-1
     signal_u1{1,i+1}=strcat('In2_u',num2str(i));    
end
signal = [signal_u0,signal_u1];

%% Initializing CBS object parameters

% Input ranges
fprintf('\n Range of throttle is [35,100]\n')
fprintf('Range of break is [0,40] \n')
fprintf('\n Grid discretization unit for both signal ranges is 4 units\n')
R1 = [35,100];
R2 = [0,40];
CBS = BrSys.copy;
CBS.SetParamRanges(signal,[ones(N1,1)*R1;ones(N2,1)*R2]);
CBS.SetEpsGridsize([4*ones(N1,1);4*ones(N2,1)]);
CBS.SetDeltaGridsize(2*CBS.epsgridsize);

fprintf('\n Press any key to continue\n ')
pause


%% Specifying STL formula
fprintf('\n The STL formula is\n ')
        f1 = STL_Formula('f1','alw(RPM[t]<2520)');
        f2 = STL_Formula('f2','ev_[0,10](speed[t]>50)');
        phi = STL_Formula('phi1','not(f1 and f2)');
        
        % Can rewrite formula: [f1 => phitest2] where
        phitest1 = STL_Formula('phitest1','ev(RPM[t]>2520)'); % not(f1)
        phitest2 = STL_Formula('phitest2','alw_[0,10](speed[t]<50)'); % not(f2)
disp(phi)

fprintf('\n Press any key to continue\n')
pause

%% Setting falsification method and parameters
msg1 = sprintf('\nChoose a falsification method\n');
msg2 = sprintf('Press 1 : Classification guided sampling\n');
msg3 = sprintf('Press 2:  Pseudo random sampling\n');
msg4 = sprintf('Press 3:  Global_nelder_mead\n');
msg5 = sprintf('Press 4:  CMA-ES\n');
msg6 = sprintf('Press 5:  simulannealbnd\n');

a = input([msg1,msg2,msg3,msg4,msg5,msg6]);

switch a 
    
    case 2
        time_lim = 2000; 
        fprintf('\n Time limit of computation is %d seconds\n',time_lim)
        snap_grid = 'y';
        switch snap_grid
            case 'y'
            CBS.SetSnapToGrid(true);
            case 'n'
            CBS.SetSnapToGrid(false);
            otherwise
                error('no epsilon resolution specified')
        end
        max_sim = inf;
        fprintf('\n Choose one of the following seeds for pseudorandom sampling:\n')
        r = input('0, 5000, 10000 or 15000\n');
        rng(r,'twister');  
        tic
        Out = StatFalsify(CBS,phi,w_rob,max_sim,max_sim,time_lim);
        time = toc;
        fprintf('Computation time = %f seconds \n',time);
        
    case 1 
        time_lim = inf;
        snap_grid = 'y';
        switch snap_grid
            case 'y'
            CBS.SetSnapToGrid(true);
            case 'n'
            CBS.SetSnapToGrid(false);
            otherwise
                error('no epsilon resolution specified')
        end 
        
        w_rob = 0.5;
        fprintf('\n Weightage to robustness information is %f\n',w_rob)
        
        max_sim = 1500;
        fprintf('\n Limit on number of simulations during global search is %d.\n',max_sim) 
        
        init_sim = 70;
        fprintf('\n Threshold number of samples for classification is %d\n ',init_sim)                
        
        
        fprintf('\n Press any key to continue\n')
        pause
        
        fprintf('\n Choose one of the following seeds for pseudorandom sampling:\n')
        r = input('0, 5000, 10000 or 15000\n');
        rng(r,'twister');  
        timervar_1 = tic;
        Out = StatFalsify([],CBS,phitest1,w_rob,init_sim,max_sim,time_lim);
        time_1 = toc(timervar_1);  
        fname = ['cl',num2str(r)];
        save(fname, 'AutoTrans_phitest1')
 %%     
       if isempty(Out.falsifier)
        % Resetting CBS
        clear CBS;
        CBS = BrSys.copy;
        CBS.SetParamRanges(signal_u0,ones(N1,1)*R1);
        CBS.SetParamRanges(signal_u1,ones(N2,1)*R2);
        CBS.SetEpsGridsize([4*ones(N1,1);4*ones(N2,1)]);
        CBS.SetDeltaGridsize(2*CBS.epsgridsize);

        
        % Sort region indices in ascending order of lowest robustness values
        [~,I] = sort(Out.lower_bounds.vals);
         L = Out.lower_bounds;
         X = L.pts;
         avg = mean(L.vals);
         Y = L.pts(L.vals<avg,:);
         es_time = 1000;
         timervar_2 = tic;
         
        switch (r)
            case {0, 10000}
             Z = X;
            case {5000, 15000}
             Z = Y;               
        end
            delete('var*','outcm*')
            rng(r,'twister');
            falsif_pb = FalsificationProblem(CBS, phi);  
            falsif_pb.setup_solver('cmaes');
            falsif_pb.solver_options.SaveVariables = 'off';
            falsif_pb.solver_options.Seed = r;
            falsif_pb.max_time = es_time;
            falsif_pb.x0 = Z;
            falsif_pb.solve()
            trace = falsif_pb.GetBrSet_False();
            if ~isempty(trace)
                fprintf('falsified')
                time_2 = toc(timervar_2);
                fprintf('\n Global computation time = %f seconds \n', time_1);
                fprintf('\n Local computation time = %f seconds\n', time_2);
                fprintf('\n Total computation time = %f seconds \n',time_1+time_2);
                trace.PlotSignals
            end
       end
       
      
    case 3
        %%
        time_lim = 2000; 
        fprintf('\n Time limit of computation is %d seconds\n',time_lim)
        falsif_pb = FalsificationProblem(CBS, phi);  
        falsif_pb.setup_solver('global_nelder_mead');
        falsif_pb.max_time = time_lim;
        timervar_2 = tic;
        falsif_pb.solve()
        trace = falsif_pb.GetBrSet_False();
        time = toc(timervar_2);
        fprintf('Computation time = %f seconds \n',time);
      
    case 4 
        %%
        time_lim = 2000; 
        fprintf('\n Time limit of computation is %d seconds\n',time_lim)
        delete('var*','outcm*')
        falsif_pb = FalsificationProblem(CBS, phi);  
        falsif_pb.setup_solver('cmaes');
        falsif_pb.solver_options.SaveVariables = 'off';
        fprintf('\n Choose one of the following seeds for pseudorandom sampling:\n')
        r = input('0, 5000, 10000 or 15000\n');
        falsif_pb.solver_options.Seed = r;
        falsif_pb.max_time = time_lim;
        timervar_2 = tic;
        falsif_pb.solve()
        trace = falsif_pb.GetBrSet_False();
        time = toc(timervar_2);
        fprintf('Computation time = %f seconds \n',time);
        if ~isempty(trace)
            trace.PlotSignals
        end
    case 5
        %%
        time_lim = 2000; 
        fprintf('\n Time limit of computation is %d seconds\n',time_lim)
        falsif_pb = FalsificationProblem(CBS, phi);  
        falsif_pb.setup_solver('simulannealbnd');
        falsif_pb.max_time = time_lim;
        timervar_2 = tic;
        falsif_pb.solve()
        trace = falsif_pb.GetBrSet_False();
        time = toc(timervar_2);
        fprintf('Computation time = %f seconds \n',time);
    otherwise
            error('no input')
end     
