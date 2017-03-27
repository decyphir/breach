function Sp = Straj2pts(S)
    
  X = cat(2,S.traj{1}.X);
    
  if (size(X,1)>numel(S.pts(:,1)))
    X = X(1:end-1, :);
  end
  
  Sp.dim = S.dim;
  Sp.pts= X;

  if isfield(S.traj, 'Expa') 
    Sp.epsi= cat(2,S.traj{1}.Expa);
  else
    Sp.epsi= repmat(S.epsi(:,1),1,size(X,2));
  end
  
