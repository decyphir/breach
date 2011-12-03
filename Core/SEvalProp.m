function [S,val] =  SEvalProp(Sys,S,props, tau, ipts)
%
%   SEVALPROP Eval property for previously computed trajectories
%  
%   Usage: [Pf val] = SEvalProp(Sys, Ptraj ,prop,tau, ipt)
%   
%   Inputs: 
%   
%    - Sys      system    
%    - Ptraj    param set with trajectories
%    - prop     property(ies)  
%    - tau      time instant(s) when to estimate properties
%    - ipts     trajectories for which to eval properties 
% 
%   Outputs: 
%  
%    - Pf       param set with prop_values field 
%    - val      quantitative satisfaction of properties
%   
         
  if (~exist('ipts')||isempty(ipts))
    ipts = 1:numel(S.traj);
  end

  if ~isfield(S,'props')
    S.props = [];
    npb =0;
  else
    npb = numel(S.props);
  end
  
  if ~isfield(S,'props_names')
    S.props_names = {} ;		
  end  
  

 %if (~exist('props')||isempty(props))
 %   props = S.props;
 % else    
 %   npb = numel(S.props);
 %   S.props = [S.props props];
 % end
 
  if (~exist('tau')||isempty(tau))
    tau0=[];
  else
    tau0 = tau;
  end
     
  for np = npb+1:numel(props)+npb
    prop = props(np-npb);
    prop_name =  get_id(prop);
    iprop = find_prop(S,prop_name);

    if ~iprop      
      S.props_names= {S.props_names{:} get_id(prop)};
      S.props= [S.props prop];
      iprop = numel(S.props_names);      
    end    

    prop = QMITL_OptimizePredicates(Sys,prop);
    fprintf(['Checking ' prop_name  '\n[             25%%           50%%            75%%               ]\n ']);
    iprog =0;
    for i = ipts
      while (floor(60*i/numel(ipts))>iprog)
        fprintf('^');
        iprog = iprog+1;
      end
      
      traj = S.traj(i);
      if (~isempty(tau0))        
        S.props_values(iprop,i).tau = tau0;
        S.props_values(iprop,i).val = QMITL_Eval(Sys,prop,traj, tau0);
        val(i) =  S.props_values(iprop,i).val(1);
      else
        tau = traj.time; 
        S.props_values(iprop,i).tau = traj.time;
        S.props_values(iprop,i).val = QMITL_Eval(Sys,prop, traj, tau);         
        val(i) =  S.props_values(iprop,i).val(1);
      end      
    end  
  end
     
  fprintf('\n');
  
function i = find_prop(S,st)

  i=0;
  for k = 1:numel(S.props_names)
    if strcmp(st,S.props_names{k})
      i = k;
      return;
    end    
  end   