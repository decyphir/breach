function  phi = RecoverFormula(phi_old)
% RECOVERFORMULA Tries to convert an STL formula created with an older version of Breach
%
% Synopsis:  phi = RecoverFormula(phi_old) 
%
% where phi is a valid QMITL formula and phi_old is a struct converted from
% a QMITL formula created with an older version of Breach.
%

if strcmp(phi_old.type, 'predicate')
    st = phi_old.st;
else
    st = form_string(phi_old);
end
    phi = QMITL_Formula(phi_old.id, st);
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
        st = cell(1,numel(phi.phin));
        for ii=1:numel(phi.phin)
            st{ii} = form_string(phi.phin(ii));
        end
        st = sprintf('%s, ',st{:});
        st = ['andn(',st(1:end-2),')'];
        
    case 'always'
        st = form_string(phi.phi);
        
        intst = phi.interval;
        try
            I = eval(phi.interval);
        catch
            I = [0 0];
        end
        
        if(I==[0 inf])
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
        
        if(I==[0 inf])
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
        
        
        if(I==[0 inf])
            st = ['(' st1 ') until (' st2 ')'];
        else
            st = ['(' st1 ')' ' until_' intst ' (' st2 ')'];
        end
        
        
end

end

function intst = form_interval_string(interval)

intst = [ '[' num2str(interval(1)) ', ' num2str(interval(2)) ']' ];

end

