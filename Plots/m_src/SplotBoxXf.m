function SplotBoxXf(S,c,alpha)
  
  nb_pts = numel(S.pts(1,:));

  if (nb_pts==0)
    return;
  end
    
  if (nargin == 1)
    c = 'g';
    alpha = .1;
    depth =0;
  end
   
  depth=1;
  Sp = S;
  
  while(isfield(Sp,'child'))
    Sp = Sp.child;
    depth = depth+1;
  end

  Sp = S;
  hold on;
 
  if (isfield(S,'plot_proj'))
    proj = S.plot_proj;
  else
    switch (size(S.Xf,1))
     case {2}
      proj=[ 1 2 ];
     otherwise
      proj = [1 2 3];
    end
  end
    
  switch(numel(proj))
   
  case 2
       
    for i=1:depth      
      
      nb_pts= size(Sp.pts,2);

      for k=1:nb_pts
        DX = Sp.traj{k}.Expa(proj,end)';
        X = Sp.traj{k}.X(proj,end)'-DX;        
        rect(X,2*DX,c,alpha);
      end
      
      if (i<depth)
        Sp = Sp.child;
      end
      
    end
    
   otherwise

    for i=1:depth      
      
      nb_pts= size(Sp.pts,2);

      for k=1:nb_pts
        DX = Sp.traj{k}.Expa(proj,end)';
        X = Sp.traj{k}.X(proj,end)'-DX;        
        voxel(X,2*DX,c,alpha);
      end
      
      if (i<depth)
        Sp = Sp.child;
      end     
    
    end
    
  end
  
