function out_signal_names = get_out_signal_names(phi)
%GET_OUT_SIGNAL_NAMES gets all output signals in a formula
%
%  syntax : out_signal_names = get_out_signal_names(phi)
%
%  where out_signal_names is an array of strings of the form 
%  sig_name_1, ..., sig_name_n where each sig_name_i 
%  is the name of an output signal

    out_signal_names = phi.out_signal_names;
end
