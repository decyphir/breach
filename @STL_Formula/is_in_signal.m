function verdict = is_in_signal(phi, signal_name)
%IS_IN_SIGNAL check whether a signal (given by its name) is declared as
%  an input
%
%  syntax : verdict = is_in_signal(phi, signal_name)
%
%  where signal_name is a string containing the name of the signal.
%  returns 1 if signal_name is declared as an input signal and 0 
%  otherwise.

  verdict = 0;
  for(i=1:length(phi.in_signal_names))
     name = phi.in_signal_names{i};
     if (strcmp(name, signal_name) == 1)
         verdict = 1;
         break;
     end
  end
end
