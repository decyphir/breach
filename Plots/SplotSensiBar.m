function [h Mend opts] =  SplotSensiBar(Sys, S, ipts, opts)
% SPLOTSENSIBAR plots an histogram of sensitivities
%
%  Synopsis:  [h Mend opts] = SplotSensiBar(Sys, S, ipts [, opts])
%
%  Plots 3d histogram of logarithmic sensitivities of state variables iX
%  w.r.t. parameters iP at a given time t for a parameter vector ipts in
%  S. iX, iP and t are provided through an input dialog box, except
%  when args is given
%           
%  opts has the following fields : args, props and taus. 
%    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check inputs 
 
  if (isempty(S.pts))
    disp('S empty !');
    return
  end

  if (~exist('ipts'))
    ipts=[];    
  end
  
  if (isempty(ipts))
    ipts = 1:size(S.pts,2);
  end
  
  % deal with options

  % default values
  open_dialog = 1;    
  cutoff =0;

  if (exist('opts'))
    
    % do we need a dialog box to enter arguments (default: yes)
    if isfield(opts,'open_dialog')
      open_dialog = opts.open_dialog;
    end
   
    % what type of computation do we do (default: average)
    
    if isfield(opts,'stat_type')
      stat_type = opts.stat_type;
    end
  
    % cut off limit in percentage of the highest value (default: 0)
    
    if isfield(opts,'cutoff')
      cutoff = opts.cutoff;
    end
    
  else
    opts=[];
  end
    
  if (open_dialog)
  
    % deal with dialog box for histogram axes 
    try 
      [iX iP tspan args] = GetArgSensiBar(Sys.DimX, Sys.ParamList, opts.args);
    catch
      try
        [iX iP tspan args] = GetArgSensiBar(Sys.DimX, Sys.ParamList);
      catch
        s = lasterror;
        warndlg(['Problem: ' s.message ', giving up'] );
        close;
        return;
      end
    end
  else
    iX = opts.args.iX;
    iP = opts.args.iP;
    tspan = opts.args.tspan;
  end
      
  % for properties evaluation 
    
  if (isfield(opts, 'props'))
    props = opts.props;
  else
    props={};
  end
  
  if isfield(opts, 'taus')
    taus = opts.taus;
  else
    taus = zeros(numel(props),1);
  end
  
  if (numel(taus)==1)
    taus = repmat(taus,numel(props),1);
  end
  
  if (~isnumeric(iX))
    NiX = iX;
    iX = [];
    for i = 1:numel(NiX)
      ind = FindParam(S,NiX{i});
      iX(i) = ind;
    end    
  end
  
  if (~isnumeric(iP))
    NiP = iP;
    iP = [];
    for i = 1:numel(NiP)
      ind = FindParam(S,NiP{i});
      iP(i) = ind;
    end      
  end
  
  % From now on I shoud have Sys, ipts, tspan, iX, iP, prop, and taus 
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   Recompute trajectories if needed
 
  if (numel(tspan)==1)
    stat_type = 'aver_end';  % compute sensitivities at a fixed time
    tspan = {0 tspan};
  else
    stat_type = 'aver_max'; 
  end
 
  Ptmp = CreateSampling(Sys, iP);
  Ptmp.pts = S.pts(:,ipts);
  S = ComputeTrajSensi(Sys, Ptmp,  tspan);   
  
  %  Compute the values for the histo bars 
  
  Mend = zeros(numel(iX)+numel(props), numel(iP));
  
  switch (stat_type)
    case {'aver_sum'},  % integrate local sensitivities over the trajectories and average
    
    % Compute bars for variable sensitivities
            
    for i = 1:numel(ipts)
      
      traj = S.traj(i);
      time = traj.time;       
      
      for j = 1:numel(iX)        
        
        for k = 1:numel(iP)
          
          is = (find(S.dim==iP(k))-1)*size(traj.X,1)+iX(j);      

          dx = traj.XS(is,:);  % dX/dp[t]          
          x = traj.X(iX(j),:);  % X[t]        
          
          % replace zeros by small quantities
          ind = find(abs(x)<1e-16);        
          x(ind) = sign(x(ind))*1e-16;       
          x(x==0) = 1e-16;
          
          p = traj.param(iP(k));    % p
          xs = (dx*p)./abs(x);

          XS =  trapz(time, xs); % computes the sum
          
          % Compute the sum  
          Mend(j,k) = Mend(j,k)+XS;      
        end        
      end
    end  % end i = ipts 

    Mend = Mend/numel(ipts);
    
    % Compute bars for properties sensitivities
        
   case {'aver_end'}, 
    
    % Compute bars for variable sensitivities
            
    for i = 1:numel(ipts)
      
      traj = S.traj(i);
      time = traj.time;       
      
      for j = 1:numel(iX)        
        
        for k = 1:numel(iP)
          
          is = (find(S.dim==iP(k))-1)*size(traj.X,1)+iX(j);      

          dx = traj.XS(is,end);  % dX/dp[t]          
          x = traj.X(iX(j),end);  % X[t]        
          
          % replace zeros by small quantities
          ind = find(abs(x)<1e-16);        
          x(ind) = sign(x(ind))*1e-16;       
          x(x==0) = 1e-16;
          
          p = traj.param(iP(k));    % p
          xs = (dx*p)./abs(x);
          
          % Compute the average 
          Mend(j,k) = Mend(j,k)+xs;      
        end        
      end
    end  % end i = ipts 

    Mend = Mend/numel(ipts);
    
    % Compute bars for properties sensitivities
        
    for j = 1:numel(props)        
        
      for k = 1:numel(iP)
        
        [p x dx] = QMITL_SEvalDiff(Sys, props{j}, S,  tspan, iP(k), taus(j));
                    
        % replace zeros by small quantities
        ind = find(abs(x)<1e-16);        
        x(ind) = sign(x(ind))*1e-16;       
        x(x==0) = 1e-16;
                       
        xs = (dx.*p)./abs(x);
        
        % Compute the average 
        Mend(numel(iX)+j,k) = mean(xs);      
      
      end        
    
    end                   
    % end case aver_end
   
   case {'aver_max'}, % 
    
    % Compute bars for variable sensitivities
            
    for i = 1:numel(ipts)
      
      traj = S.traj(i);
      time = traj.time;       
      
      for j = 1:numel(iX)        
        
        for k = 1:numel(iP)
          
          is = (find(S.dim==iP(k))-1)*size(traj.X,1)+iX(j);      
          %[dx idx] = max(abs(traj.XS(is,:))); 
                             
          dx = traj.XS(is, :);  % dX/dp[t]          
          x = traj.X(iX(j),:);  % X[t]
          
          % replace zeros by small quantities
          ind = find(abs(x)<1e-16);        
          x(ind) = sign(x(ind))*1e-16;       
          x(x==0) = 1e-16;
          
          p = traj.param(iP(k));    % p
          xs = (dx*p)./abs(x);
          [dx idx] = max(abs(xs));
          xs = xs(idx);
          
          % Compute the average 
          Mend(j,k) = Mend(j,k)+xs;      
        end        
      end
    end  % end i = ipts 

    Mend = Mend/numel(ipts);
    
    % Compute bars for properties sensitivities
        
    for j = 1:numel(props)        
        
      for k = 1:numel(iP)
        
        [p x dx] = QMITL_SEvalDiff(Sys, props{j}, S,  tspan, iP(k), taus(j));
                    
        % replace zeros by small quantities
        ind = find(abs(x)<1e-16);        
        x(ind) = sign(x(ind))*1e-16;       
        x(x==0) = 1e-16;
                       
        xs = (dx.*p)./abs(x);
        
        % Compute the average 
        Mend(numel(iX)+j,k) = mean(xs);      
      
      end        
    
    end                   
    
  end % end switch
  
  % Cut off negligible values
  M = max(max(abs(Mend)));
  Mend(abs(Mend)<cutoff*M) = 0;
  
  h =plot_histo(Mend,S, iX,props,iP);
  
function h = plot_histo(Mend,S,iX, props, iP)
  
  h = figure;    
  h = bar3(Mend,0.5,'detached');
  
  % Ticks labels
 
  xtick_labels = {};
  ytick_labels = {};
  
  for j = 1:numel(iX)
    ytick_labels = {ytick_labels{:}  S.ParamList{iX(j)}};          
  end
  
  for j = 1:numel(props)
    ytick_labels = {ytick_labels{:}  get_id(props{j})};          
  end
  
  for k = 1:numel(iP)

    xlabel = S.ParamList{iP(k)};
    if (iP(k)<= S.DimX)
      xlabel = [xlabel '(0)'];
    end
    xtick_labels = {xtick_labels{:}, xlabel };                    

  end
        
  set(gca, 'XTick', 1:size(Mend,2));
  set(gca, 'YTick', 1:size(Mend,1));
  set(gca, 'XTickLabel',  xtick_labels );
  set(gca, 'YTickLabel',  ytick_labels );
  axis([0 size(Mend,2)+1 0 size(Mend,1)+1]);
  
  hx = get(gca, 'xlabel')
  set(hx, 'Interpreter','none');
  hy = get(gca, 'ylabel')
  set(hy, 'Interpreter','none');
  
  shading interp;
  colormap cool;
  colorbar;
  for i = 1:length(h)
    zdata = get(h(i),'Zdata');
    set(h(i),'Cdata',zdata);
    set(h,'EdgeColor','k');
  end