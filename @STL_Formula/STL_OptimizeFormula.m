function phi = STL_OptimizeFormula(phi)
%STL_OptimizeFormula optimizes recursivly a formula.
% 
% Synopsis: phi = STL_OptimizeFormula(phi)
% 
% Input:
%  - phi : the formula to optimize
% 
% Output:
%  - phi : the optimized formula
% 
% Example:
%
% 
%See also
%

if strcmp(phi.type,'andn')
    phin = [];
    updated = false; % true if a formula in phin is type "and" or "andn"
    for ii=1:numel(phi.phin) % for each subformula of phi
        if strcmp(phi.phin(ii).type,'and') % if it is a "and", we add its two subformulas to phin
            phin = [phin, phi.phin(ii).phi1, phi.phin(ii).phi2]; %#ok<AGROW>
            updated = true;
        elseif strcmp(phi.phin(ii).type,'andn') % if it is a "andn", we add all its subformula to phin
            phin = [phin, phi.phin(ii).phin]; %#ok<AGROW>
            updated = true;
        else
            phin = [phin, phi.phin(ii)]; %#ok<AGROW>
        end
    end
    phi.phin = phin;
    if(updated) % if some changes are done, we apply the function again (if some subformula in phin are "and" or "andn")
        phi = STL_OptimizeFormula(phi);
    else
        for ii=1:numel(phi.phin)
            phi.phin(ii) = STL_OptimizeFormula(phi.phin(ii)); % otherwize, we go down in the tree
        end
    end
    
elseif strcmp(phi.type,'and') % we concatenate only one subformula (if appliable)
    if strcmp(phi.phi1.type,'and')
        phi.type = 'andn'; % change the subformula type
        phi.phin = [phi.phi1.phi1, phi.phi1.phi2, phi.phi2]; % 
        phi.phi1 = []; % clear unused field (what about remove them?)
        phi.phi2 = [];
        phi = STL_OptimizeFormula(phi); % apply the function again (if some subformula in phin are "and" or "andn")
        
    elseif strcmp(phi.phi1.type,'andn')
        phi.type = 'andn';
        phi.phin = [phi.phi1.phin, phi.phi2];
        phi.phi1 = [];
        phi.phi2 = [];
        phi = STL_OptimizeFormula(phi);
        
    elseif strcmp(phi.phi2.type,'and')
        phi.type = 'andn';
        phi.phin = [phi.phi1, phi.phi2.phi1, phi.phi2.phi2];
        phi.phi1 = [];
        phi.phi2 = [];
        phi = STL_OptimizeFormula(phi);
        
    elseif strcmp(phi.phi2.type,'andn')
        phi.type = 'andn';
        phi.phin = [phi.phi1, phi.phi2.phin];
        phi.phi1 = [];
        phi.phi2 = [];
        phi = STL_OptimizeFormula(phi);
        
    else % no reduction is possible, so optimize subformula
        phi.phi1 = STL_OptimizeFormula(phi.phi1);
        phi.phi2 = STL_OptimizeFormula(phi.phi2);
        
    end
    
else % neither an "and" or an "andn", so optimize subformula
    if ~isempty(phi.phi)
        phi.phi = STL_OptimizeFormula(phi.phi);
    elseif ~isempty(phi.phi1)
        phi.phi1 = STL_OptimizeFormula(phi.phi1);
    elseif ~isempty(phi.phi2)
        phi.phi2 = STL_OptimizeFormula(phi.phi2);
    end
    
end

end
