function BreachReport(s)
  % BreachReport log journal for Breach events, beginning with errors   
    
  global BreachGlobOpt
  
  if (~isfield(BreachGlobOpt,'log'))
    BreachGlobOpt.log.errors= {};
  else
    if (~isfield(BreachGlobOpt.log,'errors'))
      BreachGlobOpt.log.errors= {};
    end
  end
  
  BreachGlobOpt.log.errors = {BreachGlobOpt.log.errors{:}, s};
 
  if numel(BreachGlobOpt.log.errors)>50
    BreachGlobOpt.log.errors = {BreachGlobOpt.log.errors{end-49:end}};
  end
  
    
 
