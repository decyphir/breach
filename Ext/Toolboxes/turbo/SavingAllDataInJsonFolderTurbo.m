function SavingAllDataInJsonFolderTurbo(lower_bound, upper_bound, maxNumberOfIterations)

%% Initialize the scenario: Input parameterization for falsification problem in a "json" file
inputsRange = [lower_bound, upper_bound];
numberOfInputs = length(inputsRange(:,1));

%% Creating the scenario
scenario = {};
scenario.optimization_iterations = maxNumberOfIterations;
scenario.number_dimensions = numberOfInputs;

scenario.input_parameters = {};

%% Generating input parameters: x1,x2,...,xn and filling their values 
for i = 1: numberOfInputs
   eval(['x' num2str(i) ' = '  mat2str(inputsRange(i,:)) ';'])
end

%% Adding input parametrs to main scenario
for i =  1:numberOfInputs
   eval(['scenario.input_parameters.x' num2str(i) ' = x' num2str(i) ';'])
end

%% Creating "json" file and save all assigned values
jsonFileName = sprintf('scenario.json');
openJSONFile = fopen(jsonFileName,'w');  
encodedJSON = jsonencode(scenario); 

%% Cleaning the "json" file to make it more readable and every data to be in one line
encodedJSONCommaClean = strrep(encodedJSON, ',"', sprintf(',\r"'));
encodedJSONOpenBracketClean = strrep(encodedJSONCommaClean, '{"', sprintf('{\r"'));
encodedJSONCloseBracketClean = strrep(encodedJSONOpenBracketClean, '}', sprintf('\r}'));
encodedJSONColonBracketClean = strrep(encodedJSONCloseBracketClean, ':{', sprintf(':\r{'));
encodedJSONColonBracketCleanBracket_begin = strrep(encodedJSONColonBracketClean, '["', sprintf('"'));
encodedJSONColonBracketCleanBracket_end = strrep(encodedJSONColonBracketCleanBracket_begin, '"]', sprintf('"'));

fprintf(openJSONFile, encodedJSONColonBracketCleanBracket_end); 
fclose(openJSONFile);


end