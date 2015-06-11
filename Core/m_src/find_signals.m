function [log_sig] = find_signals(mdl) 
%  
% Find output signals and signals with a name in a Simulink model  
%  
        
  % find logged signals
  l = find_system(mdl,'FindAll','on', 'type','line');

  log_sig= {};
  for i = 1:numel(l)
    nm = get(l(i),'Name');
    if (~isempty(nm))     
      if (get(l(i),'DataLogging'))        
        log_sig = {log_sig{:}, nm};       
      end      
    end
  end
         
  log_sig = unique(log_sig);


