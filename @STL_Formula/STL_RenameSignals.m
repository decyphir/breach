function phi_new = STL_RenameSignals(phi,varargin)
%
%  Syntax: phi = STL_RenameSignals(phi, 'x', 'signal1', 'y', 'signal2',...)
%
%  In phi, rename 'x[t]' with 'signal1[t]' and 'y[t]' into 'signal2[t]', etc
%

global BreachGlobOpt

phi_new = phi;
if isequal(phi.type, 'predicate')
    
    st_phi = disp(phi,0);
    def_p=get_params(phi);
    phi_id = get_id(phi);
    phi_new_id= phi_id;
    
    
    while numel(varargin)>=2
        sig_to_rep= varargin{1};
        sig_to_rep_with = varargin{2};
        if numel(varargin)>2
            varargin = varargin(3:end);
        else
            varargin = {};
        end
        
        old = ['\<' sig_to_rep '[' ];
        new =[sig_to_rep_with '[' ];
        st_phi = regexprep(st_phi, old,new);
        phi_new_id = [phi_new_id '_' sig_to_rep_with];
        
    end
    %    phi_new_id = MakeUniqueID(phi_new_id, BreachGlobOpt.STLDB.keys);
    phi_new = STL_Formula(phi_new_id,st_phi);
    phi_new = set_params(phi_new, def_p);
    
else
    if ~isempty(phi.phi)
        phi_new.phi= STL_RenameSignals(phi.phi, varargin{:});
    end
    if ~isempty(phi.phi1)
        phi_new.phi1= STL_RenameSignals(phi.phi1, varargin{:});
    end
    if ~isempty(phi.phi2)
        phi_new.phi2= STL_RenameSignals(phi.phi2, varargin{:});
    end   
end

end
