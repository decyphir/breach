function [S,val] =  SplotSat(Sys,S,props, depth, tau, ipts)
%
%   SplotSat Plots the satisfaction function wrt time of properties for computed trajectories
%  
%   Usage: [Pf val] = SplotSat(Sys, Ptraj ,prop [, depth, tau,  ipt])
%   
%   Inputs: 
%   
%    - Sys        system    
%    - Ptraj      param set with trajectories
%    - prop       property(ies)  
%    - depth      computes and plot satisfaction of subformulas up to depth
%                 (default: 0)
%    - tau        time instant(s) when to estimate properties
%    - ipts       trajectories for which to eval properties 

%  
%   Outputs: 
%  
%    - Pf       param set with prop_values field 
%    - val      quantitative satisfaction of properties
%   
         
% check arguments 

      
  if (~exist('ipts')||isempty(ipts))
    ipts = 1:size(S.pts,2);
  end

  if (~exist('depth'))
    depth = 0;
  end
 
  if (~exist('tau')||isempty(tau))
    tau=[];
  end
 
 if isstr(props)
    QMITL_Formula('phi_tmp__', props );
    [S,val] =  SplotSat(Sys,S, phi_tmp__ , depth, tau, ipts)
    evalin('base', 'clear phi_tmp__'); 
    return;
 end

  if (depth>0)
    nprops = [];
    for i = 1:numel(props)
      nprops =   [nprops QMITL_Break(props(i), depth) ];      
    end    
    props = nprops;
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

  if ~isfield(S,'traj_ref')
    S.traj_ref =1:numel(S.traj);
  end
  
  
  % do things
    
  %% setup plots if needed
  
  nb_prop = numel(props);
  if (isfield(Sys,'time_mult'))
      time_mult = Sys.time_mult;
  else
      time_mult=1;
  end
  
  for np = npb+1:numel(props)+npb
     
    prop = props(np-npb);
    prop_name =  get_id(prop);
    iprop = find_prop(S,prop_name);
    
    subplot(nb_prop, 1, np-npb);     
    hold on;
   
    title(disp(prop), 'Interpreter','none');
   
    
    if ~iprop      
      S.props_names= {S.props_names{:} get_id(prop)};
      S.props= [S.props prop];
      iprop = numel(S.props_names);      
    end    

    prop = QMITL_OptimizePredicates(Sys,prop);
    fprintf(['Checking ' prop_name  '\n[             25%%           50%%            75%%               ]\n ']);
    iprog =0;
    
    Ptmp = Sselect(S,1);
    
    for i = ipts
      while (floor(60*i/numel(ipts))>iprog)
        fprintf('^');
        iprog = iprog+1;
      end
      
      traj = S.traj(S.traj_ref(i));
      Ptmp.pts = S.pts(:,i);
      if (~isempty(tau))        
        S.props_values(iprop,i).tau = tau;
        S.props_values(iprop,i).val = QMITL_Eval(Sys,prop,Ptmp,traj, tau);
        val(np-npb,i) =  S.props_values(iprop,i).val(1);
      else
        tau = traj.time; 
        S.props_values(iprop,i).tau = traj.time;
        S.props_values(iprop,i).val = QMITL_Eval(Sys,prop,Ptmp, traj, tau);         
        val(np-npb,i) =  S.props_values(iprop,i).val(1);
      end 

      % plot property values
      
      phi_tspan = S.props_values(iprop,i).tau;
      phi_val = S.props_values(iprop,i).val;
      plot(phi_tspan*time_mult, phi_val);
      %plot([phi_tspan(1) phi_tspan(end)]*time_mult, [0 0],'-k');
      stairs(phi_tspan*time_mult, (phi_val>0)*max(abs(phi_val))/2,'-r');
      
      legend('Quant. sat', 'Bool. sat');
      grid on;
      
    
    end  
    
    fprintf('\n');
  end
xlabel('time');    
  
function i = find_prop(S,st)

  i=0;
  for k = 1:numel(S.props_names)
    if strcmp(st,S.props_names{k})
      i = k;
      return;
    end    
  end   