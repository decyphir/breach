function SplotXU(traj,projX,projU)

  
  for k = 1:numel(traj)
    
    t = traj(k).time;
    n = numel(projX);
  
    if (nargin==2)
      projU=[];
    end
  
    nu = numel(projU);
  
    for i=1:n
      
      subplot(n+nu,1,i)
      hold on;
      plot(t,traj(k).X(projX(i),:))
      ylabel([ 'x' num2str(projX(i)) ]);

    end
  
    if (nu)
      
      for i=1:nu
      
        subplot(n+nu,1,i+n);
        hold on;
        plot(t,traj(k).U(projU(i),:));
        ylabel([ 'u' num2str(projU(i)) ]);
        
      end
      
    end
end    
