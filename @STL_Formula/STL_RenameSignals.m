function phi_new = STL_RenameSignals(phi,varargin)
%
%  Syntax: phi = STL_RenameSignals(phi, 'x', 'signal1', 'y', 'signal2',...)
%
%  In phi, rename 'x[t]' with 'signal1[t]' and 'y[t]' into 'signal2[t]', etc
%

global BreachGlobOpt


while numel(varargin)>=2
    st_phi = disp(phi);   
    sig_to_rep= varargin{1};
    sig_to_rep_with = varargin{2};
    if numel(varargin)>2
        varargin = varargin(3:end);
    else
        varargin = {};
    end
       
    old = ['\<' sig_to_rep '[' ];
    new =[sig_to_rep_with '[' ];
    phi_id = MakeUniqueID([get_id(phi) '_' sig_to_rep_with], BreachGlobOpt.STLDB.keys);
    st_new = regexprep(st_phi, old,new);
    phi_new = STL_Formula(phi_id,st_new);
    def_p=get_params(phi);
    phi_new = set_params(phi_new, def_p);
    phi = phi_new;
end
end
