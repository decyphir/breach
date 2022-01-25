function st = disp(phi, opt)
%DISP displays a formula
%
% Synopsis: st = disp(phi[, opt])
%
% Inputs:
%  - phi : the formula to display
%  - opt : (Optional, default=1) when printing predicate, if -1, the
%          function doesn't display the id, if 0, display predicates in
%          full, otherwise, display id and not predicates
%
% Outputs:
%  - st : the string
%

if ~exist('opt','var')
    opt=1;
end

if strcmp(phi.type, 'predicate')
    st = phi.st;
else
    st = form_string(phi,opt);
end

end

function st = form_string(phi, opt)
%FORM_STRING
%
% Synopsis: st = form_string(phi, opt)
%

if opt==2
    opt = -2;
elseif opt==-2
   st = phi.id;
   return
end

switch(phi.type)
    case 'predicate'
        if(opt == 0 || ~isempty(regexp(phi.id,'.+__$', 'once')))
            st = phi.st;
        else
            st = phi.id;
        end
        
    case 'not'
        st1 = form_string(phi.phi,opt);
        st = ['not (' st1 ')' ];
        
    case '=>'
        st1 = form_string(phi.phi1,opt);
        st2 = form_string(phi.phi2,opt);
        st = ['(' st1 ') => (' st2 ')'];
        
    case 'or'
        st1 = form_string(phi.phi1,opt);
        st2 = form_string(phi.phi2,opt);
        st  = ['(' st1 ') or (' st2 ')'];
        
    case 'and'
        st1 = form_string(phi.phi1,opt);
        st2 = form_string(phi.phi2,opt);
        st = ['(' st1 ') and (' st2 ')'];
        
    case 'andn'
        st = cell(1,numel(phi.phin));
        for ii=1:numel(phi.phin)
            st{ii} = form_string(phi.phin(ii),opt);
        end
        st = sprintf('%s, ',st{:});
        st = ['andn(',st(1:end-2),')'];
        
    case 'always'
        st = form_string(phi.phi,opt);
        
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
        
    case 'once'
        st = form_string(phi.phi,opt);
        
        intst = phi.interval;
        try
            I = eval(phi.interval);
        catch
            I = [0 0];
        end
        
        if(I==[0 inf])
            st = ['once (' st ')'];
        else
            st = ['once_' intst ' (' st ')'];
        end
        
    case 'historically'
        st = form_string(phi.phi,opt);
        
        intst = phi.interval;
        try
            I = eval(phi.interval);
        catch
            I = [0 0];
        end
        
        if(I==[0 inf])
            st = ['hist (' st ')'];
        else
            st = ['hist_' intst ' (' st ')'];
        end
        
    case 'eventually'
        st = form_string(phi.phi,opt);
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
        
    case 'av_eventually'
        st = form_string(phi.phi,opt);
        intst = phi.interval;
        try
            I = eval(phi.interval);
        catch
            I = [0 0];
        end
        
        if(I==[0 inf])
            st = ['av_ev (' st ')'];
        else
            st = ['av_ev_' intst ' (' st ')'];
        end
        
    case 'until'
        st1 = form_string(phi.phi1,opt);
        st2 = form_string(phi.phi2,opt);
        
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
