function z3Info = STLtoZ3(phi, z3Info)

if nargin == 1
    % No z3Info supplied
    z3Info = struct();
    z3Info.z3String = [];
    z3Info.signals = {};
    z3Info.numbersInFormula = [];
end

switch (phi.type)
    
    case 'predicate'
        phiString = disp(phi);
        
        % To format the predicate for Z3, we need to 
        % - Remove [t], [t - 10*dt] etc
        
        % Remove [t], [t - 10*dt] etc
        leftIndexes = strfind(phiString, '[');
        rightIndexes = strfind(phiString, ']');
        while ~isempty(leftIndexes)
            phiString(leftIndexes(1):rightIndexes(1)) = '';
            leftIndexes = strfind(phiString, '[');
            rightIndexes = strfind(phiString, ']');
        end
        
        numbersInFormula = getNumbersInFormula(phi, phiString);

        z3Info.z3String = phiString;
        z3Info.numbersInFormula = numbersInFormula;
        
    case 'not'
        % STL: not(phi)
        % Z3: Not(phi)
        z3InfoPhi = STLtoZ3(phi.phi, z3Info); % z3Info for phi.phi
        z3Info.z3String = ['Not(' z3InfoPhi.z3String ')'] ;
        z3Info.numbersInFormula = z3InfoPhi.numbersInFormula;
    case 'always'
        % STL: alw_[ti, tf](phi)
        % Z3: ??
        z3InfoPhi = STLtoZ3(phi.phi, z3Info); % z3Info for phi.phi
        z3Info.z3String = z3InfoPhi.z3String;
        z3Info.numbersInFormula = z3InfoPhi.numbersInFormula;
    case 'eventually'
        % STL: ev_[ti, tf](phi)
        % Z3: ??
        z3InfoPhi = STLtoZ3(phi.phi, z3Info); % z3Info for phi.phi
        z3Info.z3String = z3InfoPhi.z3String;
        z3Info.numbersInFormula = z3InfoPhi.numbersInFormula;
    case 'and'
        % STL: phi1 and phi2
        % Z3: And(phi1, phi2)
        z3InfoPhi1 = STLtoZ3(phi.phi1, z3Info); % z3Info for phi.phi1
        z3InfoPhi2 = STLtoZ3(phi.phi2, z3Info); % z3Info for phi.phi2
        z3Info.z3String = ['And(' z3InfoPhi1.z3String ', ' z3InfoPhi2.z3String ')'];
        z3Info.numbersInFormula = union(z3InfoPhi1.numbersInFormula, ...
            z3InfoPhi2.numbersInFormula);
    case 'or'
        % STL: phi1 or phi2
        % Z3: Or(phi1, phi2)
        z3InfoPhi1 = STLtoZ3(phi.phi1, z3Info); % z3Info for phi.phi1
        z3InfoPhi2 = STLtoZ3(phi.phi2, z3Info); % z3Info for phi.phi2
        z3Info.z3String = ['Or(' z3InfoPhi1.z3String ', ' z3InfoPhi2.z3String ')'];
        z3Info.numbersInFormula = union(z3InfoPhi1.numbersInFormula, ...
            z3InfoPhi2.numbersInFormula);
    case '=>'
        % STL: phi1 => phi2
        % Z3: Implies(phi1, phi2)
        z3InfoPhi1 = STLtoZ3(phi.phi1, z3Info); % z3Info for phi.phi1
        z3InfoPhi2 = STLtoZ3(phi.phi2, z3Info); % z3Info for phi.phi2
        z3Info.z3String = ['Implies(' z3InfoPhi1.z3String ', ' z3InfoPhi2.z3String ')'];
        z3Info.numbersInFormula = union(z3InfoPhi1.numbersInFormula, ...
            z3InfoPhi2.numbersInFormula);
    case 'until'
        % STL:
        % Z3
end

z3Info.signals = STL_ExtractSignals(phi);

if nargin == 1
    % Replace parameters by their actual values
    allParams = get_params(phi);
    fieldNames = fieldnames(allParams);
    for k = 1:length(fieldNames)
        paramValue = getfield(allParams, fieldNames{k});
        z3Info.z3String = strrep(z3Info.z3String, fieldNames{k}, num2str(paramValue));
    end
end
end

function numbersInFormula = getNumbersInFormula(phi, phiString)

signalsInPhi = STL_ExtractSignals(phi);

% Remove all signal names from the string
for k = 1:length(signalsInPhi)
    phiString = strrep(phiString, signalsInPhi{k}, '');
end

% Now, the only numbers in the formula are different constants, e.g. for
% the formula 3*speed[t] > 7, we have the numbers [3,7]
numbersInFormula = regexp(phiString, '\d*', 'match');
numbersInFormula = str2double(numbersInFormula); % Convert from string cell array to vector
numbersInFormula = unique(numbersInFormula);

end
