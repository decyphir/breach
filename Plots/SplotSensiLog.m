function SplotSensiLog(S,iX,iP,ipts)
% SPLOTSENSILOG    Plots trajectories logaritmic sensibilities 
%   
%  
%   Synopsis:  SplotSensiLog(S,iX,iP,ipt)
%
%     Plots trajectories logaritmic sensibilities of state variables iX
%     w.r.t. parameters iP
%   
%   Prerequisite: S has a traj field with an XS field.    
%
%   Inputs: 
%   
%    -  S        rectangular sampling set. 
%
%    -  iX       indice of the X variables for which to plot the sensitivity (idem)
%   
%    -  iP       indice of the sensitive parameter (optional, all if absent)
%
%    -  ipts     (optional) indices of the trajectories in S to plot (all
%                if absent)
%
%  
%  SEE ALSO SPLOTSENSI  
%  
  
  
  
%   Check inputs 
     
  if (isempty(S.pts))
    disp('S empty !');
    return
  end
  
  if (~isfield(S,'traj'))
    disp('No trajectory computed for this sampling')
  end
  
  if (nargin == 1)
    iX = 1:S.DimX;
    iP =S.dim;
    ipts = 1:numel(S.traj);
  elseif (nargin == 2)
    iP = S.dim;
    ipts = 1:numel(S.traj);
  elseif (nargin == 3)
    ipts = 1:numel(S.traj);
  end

 if (~isempty(iX))
    if (~isnumeric(iX))  
      NiX = iX;
      iX = [];
      for i = 1:numel(NiX)
        ind = FindParam(S,NiX{i});
        iX(i) = ind;
      end    
    end
  else
    iX = 1:S.DimX;
  end
  iX = iX(iX<=S.DimX);

  if (isempty(iP))
    iP = S.dim;    
  elseif (~isnumeric(iP))  
    NiP = iP;
    iP = [];
    for i = 1:numel(NiP)
      ind = FindParam(S,NiP{i});
      iP(i) = ind;
    end    
  end
    
 
  if (isfield(S,'traj_plot_opt'))
    opt = S.traj_plot_opt;
  else
    opt = {'b','LineWidth',1};
  end

  colors = hsv(numel(iP));

  if (isfield(S,'plot_proj'))
    proj = S.plot_proj;
  else
    proj = 1:size(S.pts,1);
  end
       
  if (isfield(S,'traj_plot_opt'))
    opt = S.traj_plot_opt;
  end
          
  for i = ipts
    
    time = S.traj{i}.time;       

    for j = 1:numel(iX)
      if (numel(iX)>1)
        subplot(numel(iX),1,j)
      end
      hold on;
      ylabel(['Sensi(' S.ParamList{iX(j)} ')'],'Interpreter','none');

      leg = {};

      for k = 1:numel(iP)

        is = (find(S.dim==iP(k))-1)*size(S.traj{i}.X,1)+iX(j);      
        dx = S.traj{i}.XS(is,:);  % dX/dp[t]        

        x = S.traj{i}.X(iX(j),:);  % X[t]        
        % replace zeros by small quantities
        ind = find(abs(x)<1e-16);        
        x(ind) = sign(x(ind))*1e-16;       
        x(x==0) = 1e-16;

        p = S.traj{i}.param(iP(k));    % p
        
        plot(time, (dx*p)./x ,opt{:},'Color',colors(k,:)); % plot
                                                           % dX/dp[t]* p/X
        leg = {leg{:}, S.ParamList{iP(k)}};
      
      end
      legend(leg{:}, 'Interpreter', 'none');
    end
  end
  
  xlabel('time','Interpreter','none');
  
  end
  
  
  
