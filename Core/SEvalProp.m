function [S,val] =  SEvalProp(Sys,S,props,tau, ipts)
%
%   Eval property for traj in S
%  
%   SEvalProp(Sys,S,prop,tau, ipt)
%   
%   Inputs: 
%   
%    - Sys       system    
%     
%    -  S        sampling set. 
%
%    -  prop     property(ies)
%   
%    -  tau    time instant(s) when to estimate properties
%    -  ipts    trajectories for which to eval properties
% 
%   Outputs: 
%  
%    - S, val
%   
      
  
  if (~exist('ipts')||isempty(ipts))
    ipts = 1:numel(S.traj);
  end

  if ~isfield(S,'props')
    S.props = [];
  end
  
  if ~isfield(S,'props_names')
    S.props_names = {} ;		
	end  
  
  if (~exist('props')||isempty(props))
    props = S.props;
  else
    
    npb = numel(S.props);
    S.props = [S.props props];
  end
 
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
      iprop = numel(S.props_names);      
    end
        
    for i = ipts
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
     
function i = find_prop(S,st)

  i=0;
  for k = 1:numel(S.props_names)
    if strcmp(st,S.props_names{k})
      i = k;
      return;
    end    
  end  

