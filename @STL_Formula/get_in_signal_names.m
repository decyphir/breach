function in_signal_names = get_in_signal_names(phi)
%GET_IN_SIGNAL_NAMES gets all input signals in a formula
%
%  syntax : in_signal_names = get_in_signal_names(phi)
%
%  where in_signal_names is an array of strings of the form 
%  sig_name_1, ..., sig_name_n where each sig_name_i 
%  is the name of an input signal

    in_signal_names = phi.in_signal_names;
end
