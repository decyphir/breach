function res = solve_turbo(this)

%% Initialization of the parameters
fcn = 'turbo_wrapper';
numberOfIterations = 1; % Counter to count the number of iterations
maxNumberOfIterations = this.max_obj_eval ;   % limit on the number of iterations
lower_bound = this.lb;  % lower bound
upper_bound = this.ub;  % upper bound
lengthOfInputVector = length(upper_bound);
xbest = [];
fbest = [];

%% Saving the necessary data in a json file, "scenarion.json", by calling SavingAllDataInJsonFolderTurbo.m
SavingAllDataInJsonFolderTurbo(lower_bound, upper_bound, maxNumberOfIterations)

% Share files between MATLAB and turbo
nameOfInputFile = 'turboToMatlab.csv';
nameOfFileToWaitFor = 'turboFinishedWriting.dummy';
nameOfOutputFile = 'matlabToTurbo.csv'; % The name to save objective value values in a csv file
nameOfTurboWaitFile = 'matlabFinishedWriting.dummy';
nameOfMatlabInitDataFile = 'initDataFromMatlabToTurbo.csv';

% Make sure the files don't exist to begin with
if exist(nameOfInputFile, 'file')
    delete(nameOfInputFile);
end
if exist(nameOfFileToWaitFor, 'file')
    delete(nameOfFileToWaitFor);
end
if exist(nameOfOutputFile, 'file')
    delete(nameOfOutputFile);
end
if exist(nameOfTurboWaitFile, 'file')
    delete(nameOfTurboWaitFile);
end
if exist(nameOfMatlabInitDataFile, 'file')
    delete(nameOfMatlabInitDataFile);
end

% Write initial data to .csv file
M = [this.solver_options.start_sample' , ...
    this.solver_options.start_function_values'];
csvwrite(nameOfMatlabInitDataFile, M);

%% Main loop:
% Running main turbo_optimization.py file in python that runs the turbo in the background

turboPythonFileLocation = which('turbo_optimization.py');
if isempty(turboPythonFileLocation)
    error('Cannot find turbo_optimization.py. Make sure to add turbo_optimization.py to MATLAB path');
end

% Note: Using java runtime to run turbo will not open up a command window,
% which is cleaner, but we also don't access the output from turbo.
% To access output from turbo, use system() to call the python script
% instead.
if this.solver_options.use_java_runtime_call
    runtime = java.lang.Runtime.getRuntime();
    process = runtime.exec(sprintf('python %s', ...
        turboPythonFileLocation));
else
    [status, commandOut] = system(sprintf('python %s &', ...
        turboPythonFileLocation));
end

%% Main loop to calculate the objective function value in Breach
while (numberOfIterations <= maxNumberOfIterations)  % While loop over the number of iteration to search falsified point
    
    while ~exist(nameOfFileToWaitFor, 'file') % While loop to wait until that input value be generated from Turbo and saved in a csv file
        pause(1);
    end
    
    inputValue = csvread(nameOfInputFile);    % Reading input values saved in a csv file
    delete(nameOfInputFile);
    delete(nameOfFileToWaitFor);
 
    if length(inputValue) ~= lengthOfInputVector  % Make sure that the number of dimensions is true
        error('Wrong number of dimensions in input file %s', nameOfInputFile);
    end
    
    objectValue = feval(fcn, inputValue, this); % Calculating the objective function value for the current input value
    
     % Check to see that objective value is negative or not
    if objectValue < 0
        xbest  = inputValue ;
        fbest = objectValue;
        break;
    end
    
    % saving the objective value in a csv file
    fileIdOut = fopen(nameOfOutputFile, 'w'); 
    fprintf(fileIdOut,'%4.4f\n', objectValue);
    fclose(fileIdOut); % Closing the csv file
    
    % Write file to tell Turbo that MATLAB is finished writing
    fclose(fopen(nameOfTurboWaitFile, 'w'));
    
    %disp("The number of evaluation is:  " + numberOfIterations) % Displaying the number of evaluation
    numberOfIterations = numberOfIterations + 1; %Increasing the number of iterations
    
end

% Saving the point with lowest objective value 
if isempty(xbest)
    xbest = this.x_best;
    fbest = this.obj_best;
end

res = struct('bestRob',[],'bestSample',[],'nTests',[],'bestCost',[],'paramVal',[],'falsified',[],'time',[]);
res.bestSample = xbest;
res.bestRob = fbest;

% Clean up files
if exist(nameOfInputFile, 'file')
    delete(nameOfInputFile);
end
if exist(nameOfFileToWaitFor, 'file')
    delete(nameOfFileToWaitFor);
end
if exist(nameOfMatlabInitDataFile, 'file')
    delete(nameOfMatlabInitDataFile);
end

if this.solver_options.use_java_runtime_call
    % Kill python process
    process.destroy();
end

end