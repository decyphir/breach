function phi_list = STL_RenameSignals(phi,varargin)
%
%  Syntax: STL_RenameSignals('x', {'signal1', 'signal2'}, 'y', {'signal1','signal3'})
%
%  Instantiates the template replacing 'x' with 'signal1' and 'signal2' and
%  'y' with 'signal1' and 'signal3'.
%

st_phi = disp(phi);
phi_list = {};

% For now we assume that we only have one signal
sig_to_rep= varargin{1};
sig_to_rep_with_list = varargin{2};

if ischar(sig_to_rep_with_list)
    sig_to_rep_with_list =  {sig_to_rep_with_list};
end

for sig_to_rep_with = sig_to_rep_with_list
    old = ['\<' sig_to_rep '[' ];
    new =[sig_to_rep_with{1} '[' ];
    phi_id = [get_id(phi.phi) '_' sig_to_rep_with{1}];
    st_new = regexprep(st_phi, old,new);
    phi_new = STL_Formula(phi_id,st_new);
    def_p=get_params(phi);
    phi_new = set_params(phi_new, def_p);
    phi_list{end+1} = phi_new;
end

if numel(phi_list)==1
   phi_list = phi_list{1};
end

end
