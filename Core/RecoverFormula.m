function  phi = RecoverFormula(phi_old)
%RECOVERFORMULA tries to convert an STL formula created with an older
% version of Breach
%
% Synopsis:  phi = RecoverFormula(phi_old)
% 
% Input:
%  - phi_old : is a STL formula or a struct converted from a STL
%              formula created with an older version of Breach.
% 
% Output:
%  - phi : is a valid STL formula
% 
%See also STL_Formula STL_ReadFile
%

if isa(phi_old,'STL_Formula')
    phi_old = struct(phi_old);
end

if strcmp(phi_old.type, 'predicate')
    st = phi_old.st;
else
    st = form_string(phi_old);
end

phi = STL_Formula(phi_old.id, st);

end

function st = form_string(phi)
%FORM_STRING
% 
% Synopsis: st = form_string(phi)
%

switch(phi.type)
    case 'predicate'
        st = phi.st;
        
    case 'not'
        st1 = form_string(phi.phi);
        st = ['not (' st1 ')' ];
        
    case '=>'
        st1 = form_string(phi.phi1);
        st2 = form_string(phi.phi2);
        st = ['(' st1 ') => (' st2 ')'];
        
    case 'or'
        st1 = form_string(phi.phi1);
        st2 = form_string(phi.phi2);
        st  = ['(' st1 ') or (' st2 ')'];
        
    case 'and'
        st1 = form_string(phi.phi1);
        st2 = form_string(phi.phi2);
        st = ['(' st1 ') and (' st2 ')'];
        
    case 'andn'
        %Expand andn because STL_Formula cannot construct formula with it
        st = form_string(phi.phin(1));
        for ii=2:numel(phi.phin)
            st = ['(' st ') and (' form_string(phi.phin(ii)) ')']; %#ok<AGROW>
        end
        
    case 'always'
        st = form_string(phi.phi);
        
        intst = phi.interval;
        try
            I = eval(phi.interval);
        catch
            I = [0 0];
        end
        
        if isequal(I,[0 inf])
            st = ['alw (' st ')'];
        else
            st = ['alw_' intst ' (' st ')'];
        end
        
    case 'eventually'
        st = form_string(phi.phi);
        
        intst = phi.interval;
        try
            I = eval(phi.interval);
        catch
            I = [0 0];
        end
        
        if isequal(I,[0 inf])
            st = ['ev (' st ')'];
        else
            st = ['ev_' intst ' (' st ')'];
        end
        
    case 'until'
        st1 = form_string(phi.phi1);
        st2 = form_string(phi.phi2);
        
        intst = phi.interval;
        try
            I = eval(phi.interval);
        catch
            I = [0 0];
        end
        
        if isequal(I,[0 inf])
            st = ['(' st1 ') until (' st2 ')'];
        else
            st = ['(' st1 ')' ' until_' intst ' (' st2 ')'];
        end
        
end

end
