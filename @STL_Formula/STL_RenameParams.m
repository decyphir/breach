function phi_new = STL_RenameParams(phi,varargin)
%
%  Syntax: phi = STL_RenameParams(phi, 'p1', 'my_p1', 'p2', 'my_p2',...)
%
%  In phi, rename 'p1' into 'my_p1' and 'p2' into 'my_p2', etc
%

global BreachGlobOpt
    
st_phi = disp(phi,0);   
def_p=get_params(phi);
phi_id = get_id(phi);
phi_new_id= phi_id;
while numel(varargin)>=2
    par_to_rep= varargin{1};
    par_to_rep_with = varargin{2};
    if numel(varargin)>2
        varargin = varargin(3:end);
    else
        varargin = {};
    end
       
    old = ['\<' par_to_rep '\>' ];
    new =par_to_rep_with;
    st_phi = regexprep(st_phi, old,new);
    phi_new_id = [phi_new_id '_' par_to_rep_with];
end
    phi_new = STL_Formula(phi_new_id,st_phi);
    phi_new = set_params(phi_new, def_p);
end
