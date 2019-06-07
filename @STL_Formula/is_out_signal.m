function verdict = is_out_signal(phi, signal_name)
%IS_OUT_SIGNAL check whether a signal (given by its name) is declared as
%  an output
%
%  syntax : verdict = is_out_signal(phi, signal_name)
%
%  where signal_name is a string containing the name of the signal.
%  returns 1 if signal_name is declared as an output signal and 0 
%  otherwise.

  verdict = 0;
  for(i=1:length(phi.out_signal_names))
     name = phi.out_signal_names{i};
     if (strcmp(name, signal_name) == 1)
         verdict = 1;
         break;
     end
  end
end
