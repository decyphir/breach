function [fv S] = SPolygon2d(Sys,S,timespan,proj)
% -
% fv = SPolygon2D(Sys,S , timespan)
%
%  Compute a structre fv with a field vertices and a field faces ready to
%  be used with patch to represent the reachable set approximated with
%  trajectories and sensitivity matrices starting from S (no refinement).
%  
%  timespan can be : 
%  
%     [t0 tf]  : plot the reachable set for interval [t0 tf], timesteps
%     are chosen by CVodes 
%
%     [t0 t1 ... tN]  : idem but timesteps are t0,t1, etc
%
%     {t0, t1, ... tN} : plot reachable set at times t0, t1, ..., tN
%     without interpolation between time instants   
%
%  proj= [i j] where x_i and x_j are the projection axes
%  
%
   
  S = ComputeTrajSensi(Sys,S, timespan);
  
  rS = S;
  rS.epsi = S.epsi*2;
  rS = RefineEstim(rS,2);
  
  n= numel(S.dim);
  nbtraj = size(S.pts,2);
  nbetraj = 2^n;
  fv.vertices = [];
  fv.faces = [];

  if ~exist('proj')  % Not used yet ... 
    proj = [1 2];    
  end

  if iscell(timespan)

    for i = 0:size(S.pts,2)-1 % loop on trajectories
      
      X = cat(1, rS.etraj(nbetraj*i+1:nbetraj*(i+1)).X);
      lengthX = size(rS.etraj(nbetraj*i+1).X,2);
      
      for j = 1:lengthX
        vertices = reshape(X(:,j),2,[])';
        %face = convhull(vertices(:,1),vertices(:,2),  {'Qt', 'Qbb','Qc', 'QbB'})+size(fv.vertices,1);
        face = convhull(vertices(:,1),vertices(:,2),  {'Qt','Qc', 'QbB'})+size(fv.vertices,1);
        fv.vertices = [fv.vertices; vertices];
        
        try
          fv.faces = [fv.faces ; face']; % try in case face has the same
                                         % number of vertices
        catch
          nb0 = size(fv.faces,2);  % adjust otherwise 
          nb1 = numel(face);
          try 
            face = [face; face(end)*ones(nb0-nb1,1)];       
            fv.faces = [fv.faces ; face'];
          catch
            fv.faces = [fv.faces repmat(fv.faces(:,end), 1,nb1-nb0)];
            fv.faces = [fv.faces ; face'];
          end        
        end         
      end
    end
    
  else % interpolate between time points
   
    for i = 0:size(S.pts,2)-1 % loop on trajectories
      
      X = cat(1, rS.etraj(nbetraj*i+1:nbetraj*(i+1)).X);
      lengthX = size(rS.etraj(nbetraj*i+1).X,2);
      
      for j = 1:lengthX-1 % loop on timespan
        vertices1 = reshape(X(:,j),2,[])';   % vertices at time t_j    
        vertices2 = reshape(X(:,j+1),2,[])'; % vertices at time t_j+1    
        vertices = [vertices1; vertices2];
        face = convhull(vertices(:,1),vertices(:,2),{'Pp'})+size(fv.vertices-nbetraj,1);
        fv.vertices = [fv.vertices; vertices1];
        try
          fv.faces = [fv.faces ; face']; % try in case face has the same
                                         % number of vertices
        catch
          nb0 = size(fv.faces,2);  % adjust otherwise 
          nb1 = numel(face);
          try 
            face = [face; face(end)*ones(nb0-nb1,1)];       
            fv.faces = [fv.faces ; face'];
          catch
            fv.faces = [fv.faces repmat(fv.faces(:,end), 1,nb1-nb0)];
            fv.faces = [fv.faces ; face'];
          end        
        end         
      end

      fv.vertices =[fv.vertices; vertices2]; % add vertices at time t_N
      
    end
  
  
  end
